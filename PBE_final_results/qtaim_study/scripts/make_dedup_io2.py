"""Deduplicate the user's located I-O BCPs (io_all.txt) by inversion orbit.
O labels are k-th oxygen in cube order (verified). For each contact (I_a, O_k),
twin = (I_{a+4}, O_{k'}=inversion-twin of O_k); canonical orbit collapses them.
Real values from the user's table; Abramov G/V/H + 4 energies. Plain .docx/.csv.
"""
import sys, os, re, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seed_bcp import load_cube, BOHR
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.section import WD_ORIENT

BOHR3 = BOHR**3; BOHR5 = BOHR**5; HA = 2625.4996
TF = 0.3 * (3.0 * np.pi**2)**(2.0 / 3.0)
SYS = [('Y_str', 'Y-MOFst'), ('Tb_str', 'Tb-MOFst'), ('Gd_str', 'Gd-MOFst'),
       ('Y_actc', 'Y-MOFact'), ('Tb_act', 'Tb-MOFact'), ('Gd_act', 'Gd-MOFact'),
       ('Y_deep', 'Y-MOFvc'), ('Tb_deep', 'Tb-MOFvc'), ('Gd_deep', 'Gd-MOFvc')]
BASE = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign'


def shifts():
    for a in (-1, 0, 1):
        for b in (-1, 0, 1):
            for c in (-1, 0, 1): yield a, b, c


def add4(x): return x + 4 if x <= 4 else x - 4


def o_twin_map(cube):
    """k -> k' where O_{k'} is nearest oxygen to inv(O_k) = (1-f)%1."""
    origin, H, n, Z, pos, rho = load_cube(cube)
    Lat = np.array([n[i] * H[:, i] for i in range(3)]); Linv = np.linalg.inv(Lat)
    fr = pos @ Linv
    Oidx = [i for i in range(len(Z)) if Z[i] == 8]   # k-th O = Oidx[k-1]
    Lat2 = Lat
    def mic(c, a): return min(np.linalg.norm(c - pos[a] + s1*Lat2[0] + s2*Lat2[1] + s3*Lat2[2]) for s1, s2, s3 in shifts())
    tw = {}
    for k, ox in enumerate(Oidx, 1):
        inv = (1 - fr[ox]) % 1
        kk = min(range(len(Oidx)), key=lambda j: mic(inv @ Lat2, Oidx[j]))
        tw[k] = kk + 1
    return tw


def parse_io(path):
    data = {}
    sys = None
    for ln in open(path):
        s = ln.strip()
        if s.startswith('###') or not s: continue
        if '-MOF' in s and not re.match(r'^I\d', s): sys = s.split('@')[-1]; continue
        if re.match(r'^I\d+-O\d', s):
            t = s.split(); m = re.match(r'I(\d+)-O(\d+)', t[0])
            data.setdefault(sys, []).append((int(m.group(1)), int(m.group(2)), float(t[1]), float(t[2]), float(t[3])))
    return data


def gvh(rho_au, lap_au):
    g = TF * rho_au**(5/3) + lap_au/6; v = lap_au/4 - 2*g; return g, v, g + v


def main():
    data = parse_io('io_all.txt')
    recs = []
    for d, sys in SYS:
        syskey = sys.split('-MOF')[-1] if False else ('MOF' + sys.split('-MOF')[1]) if '-MOF' in sys else sys
        # sys like 'Y-MOFst'; data keyed by 'st','act','vc'? io_all uses 'Y-MOFst' etc.
        contacts = data.get(sys, [])
        cube = f'{BASE}/{d}/rho17.cube'
        if not os.path.exists(cube): cube = f'{BASE}/{d}/val.cube'
        tw = o_twin_map(cube)
        rep = {}
        for a, k, dist, rhoA3, lapA5 in contacts:
            rho_au = rhoA3 * BOHR3; lap_au = lapA5 * BOHR5
            canon = (min(a, add4(a)), min(k, tw.get(k, k)))
            if canon not in rep:
                rep[canon] = (f'I{a}-O{k}', dist, rho_au, lap_au)
        for bond, dist, rho_au, lap_au in rep.values():
            g, v, h = gvh(rho_au, lap_au); e = [-0.5*v, 0.429*g, 0.68*(-v), 0.67*g]
            recs.append((sys, bond, dist, rho_au, lap_au, g, v, h, [x*HA for x in e]))
    order = ['Y-MOFst','Tb-MOFst','Gd-MOFst','Y-MOFact','Tb-MOFact','Gd-MOFact','Y-MOFvc','Tb-MOFvc','Gd-MOFvc']
    recs.sort(key=lambda r: (order.index(r[0]), r[2]))
    from collections import Counter
    print('per system:', dict(Counter(r[0] for r in recs)), 'total', len(recs))

    import csv
    with open('BCP_IO_dedup.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['System', 'Bond', 'd(A)', 'rho(e/A3)', 'lap(e/A5)', 'rho_BCP(e/bohr3)',
                    'G(au)', 'V(au)', 'H(au)', 'E=-V/2', 'E=0.429G', 'E=0.68|V|', 'E=0.67G'])
        for sys, bond, d, rho_au, lap_au, g, v, h, e in recs:
            w.writerow([sys, bond, f'{d:.2f}', f'{rho_au/BOHR3:.3f}', f'{lap_au/BOHR5:+.3f}', f'{rho_au:.4f}',
                        f'{g:.4f}', f'{v:.4f}', f'{h:.4f}'] + [f'{x:.1f}' for x in e])

    doc = Document(); sec = doc.sections[0]; sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = sec.page_height, sec.page_width
    for m in ('left_margin', 'right_margin', 'top_margin', 'bottom_margin'): setattr(sec, m, Inches(0.3))
    doc.styles['Normal'].font.name = 'Calibri'; doc.styles['Normal'].font.size = Pt(9)
    cap = doc.add_paragraph(); cr = cap.add_run('Table S. I···O bond critical points (one representative per symmetry-equivalent chain pair; side-on polyiodide···carboxylate contacts). All BCPs located by critic2.')
    cr.font.size = Pt(10)
    hdr = ['System', 'Bond', 'd (Å)', 'ρ (e/Å³)', '∇²ρ (e/Å⁵)', 'ρ_BCP (e/bohr³)',
           'G (a.u.)', 'V (a.u.)', 'H (a.u.)', '−V/2', '0.429G', '0.68|V|', '0.67G']
    widths = [1.15, 0.7, 0.4, 0.55, 0.62, 0.82, 0.5, 0.5, 0.5, 0.45, 0.5, 0.5, 0.5]

    def setcell(cell, text, size=7):
        cell.text = ''; p = cell.paragraphs[0]; p.alignment = 1
        run = p.add_run(text); run.font.size = Pt(size); run.font.name = 'Calibri'
    t = doc.add_table(rows=1, cols=len(hdr)); t.style = 'Table Grid'; t.autofit = False
    for j, h in enumerate(hdr): setcell(t.rows[0].cells[j], h, size=7)
    for sys, bond, d, rho_au, lap_au, g, v, h, e in recs:
        vals = [sys, bond, f'{d:.2f}', f'{rho_au/BOHR3:.3f}', f'{lap_au/BOHR5:+.3f}', f'{rho_au:.4f}',
                f'{g:.4f}', f'{v:.4f}', f'{h:.4f}'] + [f'{x:.1f}' for x in e]
        cells = t.add_row().cells
        for j, val in enumerate(vals): setcell(cells[j], val, size=7)
    for j, wd in enumerate(widths):
        for row in t.rows: row.cells[j].width = Inches(wd)
    doc.save('BCP_IO_dedup.docx'); print('wrote BCP_IO_dedup.docx + .csv')


if __name__ == '__main__':
    main()
