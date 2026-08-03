"""M-I (Y/Tb/Gd ... I) BCP table — recovered by seeded Newton (critic2 auto missed them).
Columns: d, rho, lap, rho_BCP, G, V, H (Abramov). NO energy correlations.
Deduplicated by inversion (1/2,1/2,1/2) orbit. Plain .docx + .csv.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seed_bcp import load_cube, fields, Interp, BOHR
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.section import WD_ORIENT

BOHR3 = BOHR**3; BOHR5 = BOHR**5; TF = 0.3 * (3.0 * np.pi**2)**(2.0 / 3.0)
MZ = {'Y_str': 39, 'Tb_str': 65, 'Gd_str': 64, 'Y_actc': 39, 'Tb_act': 65, 'Gd_act': 64,
      'Y_deep': 39, 'Tb_deep': 65, 'Gd_deep': 64}
MSYM = {39: 'Y', 65: 'Tb', 64: 'Gd'}
SYS = [('Y_str', 'Y-MOFst'), ('Tb_str', 'Tb-MOFst'), ('Gd_str', 'Gd-MOFst'),
       ('Y_actc', 'Y-MOFact'), ('Tb_act', 'Tb-MOFact'), ('Gd_act', 'Gd-MOFact'),
       ('Y_deep', 'Y-MOFvc'), ('Tb_deep', 'Tb-MOFvc'), ('Gd_deep', 'Gd-MOFvc')]
BASE = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign'


def shifts():
    for a in (-1, 0, 1):
        for b in (-1, 0, 1):
            for c in (-1, 0, 1): yield a, b, c


def newton(seed, ip, n=100, tol=1e-4):
    x = seed.copy()
    for _ in range(n):
        g = np.array([ip('g' + a, x) for a in 'xyz'])
        if np.linalg.norm(g) < tol: break
        Hx = np.array([[ip('h' + a + b, x) for b in 'xyz'] for a in 'xyz'])
        try: dx = np.linalg.solve(Hx, g)
        except np.linalg.LinAlgError: break
        st = np.linalg.norm(dx)
        if st > 0.5: dx *= 0.5 / st
        x = x - dx
    return x, np.linalg.norm(g)


def recover(cube, mz):
    origin, H, n, Z, pos, rho = load_cube(cube)
    Lat = np.array([n[i] * H[:, i] for i in range(3)]); Linv = np.linalg.inv(Lat)
    g, Hxmat, lap = fields(rho, H)
    arrs = {'rho': rho, 'lap': lap}
    for k, a in enumerate('xyz'): arrs['g' + a] = g[k]
    for ia, a in enumerate('xyz'):
        for jb, b in enumerate('xyz'): arrs['h' + a + b] = Hxmat[ia, jb]
    ip = Interp(origin, np.linalg.inv(H), n, arrs)
    fr = pos @ Linv
    M = [i for i in range(len(Z)) if Z[i] == mz]
    Iidx = [i for i in range(len(Z)) if Z[i] == 53]; Ilab = {i: k + 1 for k, i in enumerate(Iidx)}

    def mic(a, b): return min(np.linalg.norm(pos[a] - pos[b] + s1*Lat[0] + s2*Lat[1] + s3*Lat[2]) for s1, s2, s3 in shifts()) * BOHR
    reps = {}
    for m in M:
        for i in Iidx:
            if mic(m, i) > 3.7: continue
            disp = fr[i] - fr[m]; disp -= np.round(disp)
            xr, gn = newton(pos[m] + (0.5 * disp) @ Lat, ip)
            if gn > 1e-3: continue
            ds = sorted([(min(np.linalg.norm(xr - pos[a] + s1*Lat[0] + s2*Lat[1] + s3*Lat[2]) for s1, s2, s3 in shifts()), a) for a in range(len(Z))])[:2]
            if {Z[ds[0][1]], Z[ds[1][1]]} != {mz, 53}: continue
            fpos = (xr) @ Linv
            key = min(tuple(np.round(fpos % 1, 3)), tuple(np.round((1 - fpos) % 1, 3)))
            if key in reps: continue
            rho_b = ip('rho', xr); lap_b = ip('lap', xr)
            reps[key] = (f'{MSYM[mz]}-I{Ilab[i]}', mic(m, i), rho_b, lap_b)
    return list(reps.values())


def gvh(rho_au, lap_au):
    g = TF * rho_au**(5/3) + lap_au/6; v = lap_au/4 - 2*g; return g, v, g + v


def main():
    recs = []
    for d, sys in SYS:
        cube = f'{BASE}/{d}/rho17.cube'
        if not os.path.exists(cube): cube = f'{BASE}/{d}/val.cube'
        for bond, dist, rho_au, lap_au in recover(cube, MZ[d]):
            g, v, h = gvh(rho_au, lap_au)
            recs.append((sys, bond, dist, rho_au, lap_au, g, v, h))
    order = ['Y-MOFst','Tb-MOFst','Gd-MOFst','Y-MOFact','Tb-MOFact','Gd-MOFact','Y-MOFvc','Tb-MOFvc','Gd-MOFvc']
    recs.sort(key=lambda r: (order.index(r[0]), r[2]))
    from collections import Counter
    print('per system:', dict(Counter(r[0] for r in recs)), 'total', len(recs))

    import csv
    with open('BCP_MI_dedup.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['System', 'Bond', 'd(A)', 'rho(e/A3)', 'lap(e/A5)', 'rho_BCP(e/bohr3)', 'G(au)', 'V(au)', 'H(au)'])
        for sys, bond, d, rho_au, lap_au, g, v, h in recs:
            w.writerow([sys, bond, f'{d:.2f}', f'{rho_au/BOHR3:.3f}', f'{lap_au/BOHR5:+.3f}', f'{rho_au:.4f}', f'{g:.4f}', f'{v:.4f}', f'{h:.4f}'])

    doc = Document(); sec = doc.sections[0]; sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = sec.page_height, sec.page_width
    for m in ('left_margin', 'right_margin', 'top_margin', 'bottom_margin'): setattr(sec, m, Inches(0.4))
    doc.styles['Normal'].font.name = 'Calibri'; doc.styles['Normal'].font.size = Pt(9)
    cap = doc.add_paragraph()
    cr = cap.add_run('Table S. Y/Ln···I bond critical points (seeded-Newton recovery; critic2 auto-scan did not locate these). '
                     'Ionic closed-shell contacts. G, V, H from Abramov + local virial (a.u.).')
    cr.font.size = Pt(10)
    hdr = ['System', 'Bond', 'd (Å)', 'ρ (e/Å³)', '∇²ρ (e/Å⁵)', 'ρ_BCP (e/bohr³)', 'G (a.u.)', 'V (a.u.)', 'H (a.u.)']
    widths = [1.2, 0.8, 0.5, 0.6, 0.7, 0.9, 0.6, 0.6, 0.6]

    def setcell(cell, text, size=8):
        cell.text = ''; p = cell.paragraphs[0]; p.alignment = 1
        run = p.add_run(text); run.font.size = Pt(size); run.font.name = 'Calibri'
    t = doc.add_table(rows=1, cols=len(hdr)); t.style = 'Table Grid'; t.autofit = False
    for j, h in enumerate(hdr): setcell(t.rows[0].cells[j], h, size=8)
    for sys, bond, d, rho_au, lap_au, g, v, h in recs:
        vals = [sys, bond, f'{d:.2f}', f'{rho_au/BOHR3:.3f}', f'{lap_au/BOHR5:+.3f}', f'{rho_au:.4f}', f'{g:.4f}', f'{v:.4f}', f'{h:.4f}']
        cells = t.add_row().cells
        for j, val in enumerate(vals): setcell(cells[j], val, size=8)
    for j, wd in enumerate(widths):
        for row in t.rows: row.cells[j].width = Inches(wd)
    doc.save('BCP_MI_dedup.docx'); print('wrote BCP_MI_dedup.docx + .csv')


if __name__ == '__main__':
    main()
