#!/usr/bin/env python3
"""Replace trilinear 'newton' values in the 442 CP tables with critic2-spline values
refined from the same positions (tools/critic2_new, `auto seed point`)."""
import re, csv, numpy as np
BOHR = 0.5291772109
T = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/qtaim_study_442/tables'
CAMP = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442'
CIFS = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results_442/landscape_3variants/cifs'
SYS = {'Y-MOFst(poly)': ('st_vc', 'st_poly_lS5', 'l442_st_poly_lS5_relaxed.cif',
        np.array([[8.9348,0,0],[0,9.6608,0],[-3.7625602248,0,16.5958256665]])),
       'Y-MOFvc(poly)': ('st_vc', 'vc_poly_lS3', 'l442_vc_poly_lS3_relaxed.cif',
        np.array([[8.7975,0,0],[-1.4362,8.2029,0],[-4.1207,-0.9336,17.2148]])),
       'Y-MOFact(poly)': ('act', 'act_poly_lS2', 'l442_act_poly_lS2_relaxed.cif',
        np.array([[8.9348,0,0],[0,9.6608,0],[-3.7625602248,0,16.5958256665]])),
       'Y-MOFact(mol)': ('act', 'act_mol_rd', 'l442_act_mol_rd_relaxed.cif',
        np.array([[8.9348,0,0],[0,9.6608,0],[-3.7625602248,0,16.5958256665]]))}
ROWRE = re.compile(r"\(\s*3\s*,\s*-1\s*\)\s+bond\s+([0-9.+-]+)\s+([0-9.+-]+)\s+([0-9.+-]+)\s+\d+\s+\S+\s+([0-9.Ee+-]+)\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)")

def cif_atoms(sysname, el):
    s = open(f'{CIFS}/{SYS[sysname][2]}').read()
    lab = s.index('_atom_site_fract_z')
    out = []
    for l in s[lab:].splitlines():
        q = l.split()
        if len(q) >= 5 and q[1] == el:
            out.append([float(q[2]), float(q[3]), float(q[4])])
    return np.array(out)

def dmic(fa, fb, cell):
    d = np.array(fa) - np.array(fb); d -= np.round(d)
    return np.linalg.norm(d @ cell)

for fname, dpcol, bondre in (('BCP_II_full_442.csv', 'd_pair(A)', r'I(\d+)-I(\d+)'),
                             ('BCP_IO_full_442.csv', 'd(A)', r'O(\d+)-I(\d+)')):
    rows = list(csv.DictReader(open(f'{T}/{fname}')))
    nchanged = nunref = 0
    # cache refined CPs per system
    cache = {}
    out = []
    for r in rows:
        if not str(r.get('src','')).startswith('newton') or r['System'] not in SYS:
            out.append(r); continue
        sysname, (d, tag, cifname, cellA) = r['System'], SYS[r['System']]
        if tag not in cache:
            cellB = cellA / BOHR
            cps = []
            for l in open(f'{CAMP}/{d}/refine_{tag}.out', errors='ignore'):
                m = ROWRE.search(l)
                if m:
                    cps.append((np.array([float(m.group(i)) for i in (1,2,3)]),
                                float(m.group(4)), float(m.group(5))))
            cache[tag] = (cellB, cps, cif_atoms(sysname, 'I'))
        cellB, cps, If = cache[tag]
        Of = cif_atoms(sysname, 'O')
        m = re.match(bondre, r['Bond'])
        if not m:
            out.append(r); continue
        dp = float(r[dpcol])
        if fname.startswith('BCP_II'):
            ia, ib = int(m.group(1)), int(m.group(2))
            fa, fb = If[ia-1], If[ib-1]
        else:
            ii = int(m.group(2))
            fb = If[ii-1]
            # recover the O position: the CIF O-I pair whose distance equals the reported d
            cand = [(abs(dmic(o, fb, cellA) - dp), o) for o in Of]
            cand.sort()
            fa = cand[0][1]
            if cand[0][0] > 0.01:
                out.append(r); continue
        best = None
        for fr, rho, lap in cps:
            d1 = dmic(fr, fa, cellB) * BOHR
            d2 = dmic(fr, fb, cellB) * BOHR
            score = abs((d1 + d2) - dp) + abs(d1 - d2)
            if best is None or score < best[0]:
                best = (score, rho, lap)
        if best is not None and best[0] < 0.35:
            _, rho, lap = best
            r['rho_BCP(e/bohr3)'] = f"{rho:.4f}"
            r['lap(e/A5)'] = f"{lap/BOHR**5:+.3f}"
            r['rho(e/A3)'] = f"{rho/BOHR**3:.3f}"
            r['src'] = 'critic2(seedpt)'
            nchanged += 1
        else:
            r['src'] = 'newton(trilinear)'
            nunref += 1
        out.append(r)
    hdr = list(rows[0].keys())
    with open(f'{T}/{fname}', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=hdr); w.writeheader(); w.writerows(out)
    print(f'{fname}: {nchanged} rows put on critic2-spline footing, {nunref} left trilinear')
