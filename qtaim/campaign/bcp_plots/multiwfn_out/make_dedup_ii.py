"""Deduplicated I-I BCP table from critic2 bcp.out (real located values).
Extract located I-I BCPs (2-nearest-I < 1.8 A, equidistant), deduplicate by the
I_a<->I_{a+4} twin orbit, compute G/V/H + 4 energies. Writes plain .docx + .csv.
"""
import sys, os, re, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seed_bcp import load_cube, BOHR
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.section import WD_ORIENT

ROWRE = re.compile(r"\(\s*3\s*,\s*-1\s*\)\s+bond\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+\d+\s+\S+\s+([0-9.Ee+-]+)\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)")
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


def extract(cube, bcpout):
    origin, H, n, Z, pos, rho = load_cube(cube)
    Lat = np.array([n[i] * H[:, i] for i in range(3)])
    Iidx = list(np.where(Z == 53)[0]); Ilab = {i: k + 1 for k, i in enumerate(Iidx)}

    def micc(c1, c2): return min(np.linalg.norm(c1 - c2 + s1*Lat[0] + s2*Lat[1] + s3*Lat[2]) for s1, s2, s3 in shifts())
    bcps = []
    for ln in open(bcpout):
        m = ROWRE.search(ln)
        if not m: continue
        rho_b = float(m.group(4))
        if not (0.02 < rho_b < 0.07): continue
        f = np.array([float(m.group(1)), float(m.group(2)), float(m.group(3))]); cart = f @ Lat
        lap = float(m.group(5))
        ds = sorted([(micc(cart, pos[i]), i) for i in Iidx])[:2]
        if ds[0][0]*BOHR < 1.8 and ds[1][0]*BOHR < 1.8 and ds[1][0] / ds[0][0] < 1.5:  # real I-I BCP
            a, b = sorted([Ilab[ds[0][1]], Ilab[ds[1][1]]])
            d = micc(pos[ds[0][1]], pos[ds[1][1]]) * BOHR
            bcps.append((a, b, d, rho_b, lap))
    # dedup by twin orbit I_a<->I_{a+4}
    rep = {}
    for a, b, d, rho_b, lap in bcps:
        tw = tuple(sorted((add4(a), add4(b)))); canon = min((a, b), tw)
        if canon not in rep: rep[canon] = (a, b, d, rho_b, lap)
    return list(rep.values())


def gvh(rho_au, lap_au):
    g = TF * rho_au**(5/3) + lap_au/6; v = lap_au/4 - 2*g; return g, v, g + v


def main():
    recs = []
    for d, sys in SYS:
        cube = f'{BASE}/{d}/rho17.cube'
        if not os.path.exists(cube): cube = f'{BASE}/{d}/val.cube'
        for a, b, dist, rho_au, lap_au in extract(cube, f'{BASE}/{d}/bcp.out'):
            g, v, h = gvh(rho_au, lap_au)
            e = [-0.5*v, 0.429*g, 0.68*(-v), 0.67*g]
            recs.append((sys, f'I{a}-I{b}', dist, rho_au, lap_au, g, v, h, [x*HA for x in e]))
    recs.sort(key=lambda r: (['Y-MOFst','Tb-MOFst','Gd-MOFst','Y-MOFact','Tb-MOFact','Gd-MOFact','Y-MOFvc','Tb-MOFvc','Gd-MOFvc'].index(r[0]), r[2]))
    from collections import Counter
    print('per system:', dict(Counter(r[0] for r in recs)), 'total', len(recs))

    import csv
    with open('BCP_II_dedup.csv', 'w', newline='') as f:
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
    cap = doc.add_paragraph(); cr = cap.add_run('Table S. I–I bond critical points (one representative per symmetry-equivalent pair, P1) — topological descriptors and interaction energies.')
    cr.font.size = Pt(10)

    hdr = ['System', 'Bond', 'd (Å)', 'ρ (e/Å³)', '∇²ρ (e/Å⁵)', 'ρ_BCP (e/bohr³)',
           'G (a.u.)', 'V (a.u.)', 'H (a.u.)', '−V/2', '0.429G', '0.68|V|', '0.67G']
    widths = [1.15, 0.6, 0.4, 0.55, 0.62, 0.82, 0.5, 0.5, 0.5, 0.45, 0.5, 0.5, 0.5]

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
    doc.save('BCP_II_dedup.docx'); print('wrote BCP_II_dedup.docx + .csv')


if __name__ == '__main__':
    main()
