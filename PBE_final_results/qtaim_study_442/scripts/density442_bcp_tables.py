#!/usr/bin/env python3
"""442 CP tables — the qtaim_study pipeline applied to the density_442 systems.
Canonical search: critic2 `auto seed pair dist 4.0 nr 8` on plot-17 cubes (bcp_raw_*.out).
Old-pipeline logic (as in PBE_final_results/qtaim_study): critic2 rho/lap taken as-is;
BCPs assigned to their two nearest nuclei (27 MIC images); near-core filter
(d1 < 2.6 A, d2/d1 < 2) PLUS the iodine 4d-shell filter (any I nucleus closer than
1.35 A => drop: I.pbe-n-kjpaw carries 4d10 in valence, whose shell structure generates
spurious (3,-1) points); "missing" I-I BCPs recovered by Newton from bond midpoints
(trilinear values, flagged src=newton). rho_line_min = straight-line MIC profile minimum
(information column; equals rho_BCP when the bond path is straight). Abramov G/V/H and
the 4 correlation energies as in the Gamma-era tables.
"""
import sys, re
import numpy as np
sys.path.insert(0, '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/bcp_plots/multiwfn_out')
from seed_bcp import load_cube, fields, Interp, newton, BOHR

CAMP = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442'
OUT = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/qtaim_study_442'
BOHR3, BOHR5 = BOHR**3, BOHR**5
HA = 2625.4996
TF = 0.3 * (3*np.pi**2)**(2/3)
ROWRE = re.compile(r"\(\s*3\s*,\s*-1\s*\)\s+bond\s+([0-9.+-]+)\s+([0-9.+-]+)\s+([0-9.+-]+)\s+\d+\s+\S+\s+([0-9.Ee+-]+)\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)")
SYS = [('Y-MOFst(poly)', 'st_vc', 'st_poly_lS5'),
       ('Y-MOFvc(poly)', 'st_vc', 'vc_poly_lS3'),
       ('Y-MOFact(poly)', 'act', 'act_poly_lS2'),
       ('Y-MOFact(mol)', 'act', 'act_mol_rd')]

def mic(dcart, M):
    best = 1e9
    for a in (-1,0,1):
        for b in (-1,0,1):
            for c in (-1,0,1):
                r = np.linalg.norm(dcart + (np.array([a,b,c]) @ M))
                if r < best: best = r
    return best * BOHR

def analyse(label, d, tag):
    D = f'{CAMP}/{d}'
    origin, H, dims, Z, pos, rho = load_cube(f'{D}/rho17_{tag}.cube')
    M = H.T * np.array(dims, dtype=float)[:,None]   # row i = h_i*n_i ( broadcasting bug fix)
    Minv = np.linalg.inv(M)
    gx, Hx, lap = fields(rho, H)
    arrays = {'rho': rho, 'gx': gx[0], 'gy': gx[1], 'gz': gx[2],
              'hxx': Hx[0][0], 'hxy': Hx[0][1], 'hxz': Hx[0][2],
              'hyy': Hx[1][1], 'hyz': Hx[1][2], 'hzz': Hx[2][2],
              'hyx': Hx[0][1], 'hzx': Hx[0][2], 'hzy': Hx[1][2]}
    ip = Interp(origin, np.linalg.inv(H), dims, arrays)
    Iidx = [k for k in range(len(Z)) if Z[k]==53]
    print(f'[{label}] grid {dims}', flush=True)

    def micf(fa, fb):
        df = fa - fb; df -= np.round(df)
        return np.linalg.norm(df @ M) * BOHR
    def line_min(a, b, n=400):
        df = (pos[a]-pos[b]) @ Minv; df -= np.round(df)
        fb = pos[b] @ Minv
        return min(ip('rho', (fb + t*df) @ M) for t in np.linspace(0.2, 0.8, n))

    # ---- candidate CPs ----
    cands = []
    for l in open(f'{D}/bcp_raw_{tag}.out', errors='ignore'):
        m = ROWRE.search(l)
        if not m: continue
        fr = np.array([float(m.group(i)) for i in (1,2,3)])
        cands.append([fr @ M, float(m.group(4)), float(m.group(5)), 'critic2'])
    for ni, a in enumerate(Iidx):
        for b in Iidx[ni+1:]:
            if micf(pos[a] @ Minv, pos[b] @ Minv) > 4.5: continue
            dfp = (pos[a]-pos[b]) @ Minv; dfp -= np.round(dfp)
            mid = (pos[b] @ Minv + dfp/2) @ M
            x, gn = newton(mid, ip, n=80, tol=1e-4)
            if gn > 1e-3: continue
            if any(np.linalg.norm(x-c[0]) < 0.2 for c in cands): continue
            cands.append([x, ip('rho',x), ip('hxx',x)+ip('hyy',x)+ip('hzz',x), 'newton'])
    # ---- assign + filter ----
    II, IO, MI, IC = [], [], [], []
    nshell = nother = 0
    have_c2 = set()
    for x, rb, lb, src in cands:
        ds = sorted((mic(pos[k]-x, M), k) for k in range(len(Z)))
        d1, d2, k1, k2 = ds[0][0], ds[1][0], ds[0][1], ds[1][1]
        if not (d1 < 2.6 and d2/d1 < 2): continue
        # iodine 4d-shell filter: close to I *and* one-sided (a genuine X-I bond midpoint
        # has d2/d1 ~ 1; a shell artifact hugs one iodine with the next nucleus ~2x farther)
        if (Z[k1]==53 and d1 < 1.35 and d2/d1 > 1.5) or (Z[k2]==53 and d2 < 1.35 and d2/d1 > 1.5): nshell += 1; continue
        pair = {Z[k1], Z[k2]}
        rec = [micf(pos[k1] @ Minv, pos[k2] @ Minv), rb, lb, k1, k2, d1, x, src]
        if pair == {53}:
            II.append(rec)
            if src == 'critic2': have_c2.add(frozenset((k1,k2)))
        elif pair == {53,8}: IO.append(rec)
        elif pair == {53,39}: MI.append(rec)
        elif pair == {53,6}: IC.append(rec)
        else: nother += 1


    # ---- bond-path verification (QTAIM: a bond IS its bond path) ----
    def path_ends(x0, seeds, n=600):
        ends=[]
        for k0 in seeds:
            x=x0.copy()
            d0=dpoint(x)
            u=(pos[k0]-x - np.round((pos[k0]-x)@Minv)@M)
            x=x+0.15*u/np.linalg.norm(u)      # leave the CP along the candidate direction
            for _ in range(n):
                if dpoint(x)[0][0] < 0.9: break
                g=np.array([ip('g'+a,x) for a in 'xyz']); gn=np.linalg.norm(g)
                if gn<1e-8: break
                step=min(0.06, 0.02/max(gn,1e-12))
                x=x+step*g/gn
            ends.append(dpoint(x)[0][1])
        return frozenset(ends)
    def dpoint(x):
        return sorted((np.linalg.norm((pos[k]-x-np.round((pos[k]-x)@Minv)@M))*BOHR, k) for k in range(len(Z)))

    def verify(recs, name):
        ok=[]
        for r in recs:
            x=r[6]
            seeds=[r[3],r[4]]
            e=path_ends(x,seeds)
            if e == frozenset((r[3],r[4])): ok.append(r+[1])
            else: r.append(0); ok.append(r)
        nok=sum(1 for r in ok if r[-1]==1)
        print(f'    {name}: path-verified {nok}/{len(ok)} (rows kept, flag column)')
        return ok
    II=verify(II,'II'); IO=verify(IO,'IO'); IC=verify(IC,'IC'); MI=verify(MI,'MI')
    # dump newton CP positions -> critic2 seed-point refinement input (uniform spline footing)
    seeds=[]
    for r in II+IO+IC+MI:
        if len(r)>6 and r[6] is not None and frozenset((r[3],r[4])) not in have_c2:
            seeds.append((r[3],r[4],r[6]))
    if seeds:
        pts=[]
        for a,b,x in seeds:
            f=(x@Minv)%1.0
            pts.append("point x0 %.8f %.8f %.8f"%(f[0],f[1],f[2]))
        lines=[f"crystal rho17_{tag}.cube",f"load rho17_{tag}.cube id rho"]
        for i in range(0,len(pts),4):
            lines.append("auto seed "+" ".join(pts[i:i+4]))
        open(f'{D}/refine_{tag}.cri','w').write("\n".join(lines)+"\n")
        print(f'  refine_{tag}.cri: {len(pts)} seeds')
    def gvh(rb, lb):
        g = TF*rb**(5/3) + lb/6; v = lb/4 - 2*g; return g, v, g+v
    Iorder = {k: i+1 for i, k in enumerate(Iidx)}
    rows = ["System,Bond,d_pair(A),d1(A),rho(e/A3),lap(e/A5),rho_BCP(e/bohr3),rho_line_min(e/bohr3),G(au),V(au),H(au),E=-V/2,E=0.429G,E=0.68|V|,E=0.67G,BONDED(0.02),src,path_ok"]
    agree = []
    for r in sorted(II, key=lambda r: -r[1]):
        dA, rb, lb, a, b, d1 = r[:6]
        g, v, h = gvh(rb, lb)
        lm = line_min(a, b)
        if rb >= 0.02: agree.append(abs(rb-lm)/rb*100)
        src = 'critic2' if frozenset((a,b)) in have_c2 else 'newton'
        rows.append(f"{label},I{Iorder[a]}-I{Iorder[b]},{dA:.3f},{d1:.3f},{rb/BOHR3:.3f},{lb/BOHR5:+.3f},{rb:.4f},{lm:.4f},"
                    f"{g:.4f},{v:.4f},{h:.4f},{-v/2*HA:.1f},{0.429*g*HA:.1f},{0.68*abs(v)*HA:.1f},{0.67*g*HA:.1f},{int(rb>=0.02)},{src},{r[8]}")
    open(f'{OUT}/tables/BCP_II_full_442.csv','a').write("\n".join(rows)+"\n")
    rows2 = ["System,Bond,d(A),rho(e/A3),lap(e/A5),rho_BCP(e/bohr3),G(au),V(au),H(au),src"]
    oorder = {}
    for r in sorted(IO, key=lambda r: -r[1]):
        dA, rb, lb, a, b, d1 = r[:6]
        g, v, h = gvh(rb, lb)
        oi = a if Z[a]==8 else b
        ii = a if Z[a]==53 else b
        oorder.setdefault(oi, len(oorder)+1)
        rows2.append(f"{label},O{oorder[oi]}-I{Iorder[ii]},{dA:.3f},{rb/BOHR3:.3f},{lb/BOHR5:+.3f},{rb:.4f},{g:.4f},{v:.4f},{h:.4f},{r[7]}")
    open(f'{OUT}/tables/BCP_IO_full_442.csv','a').write("\n".join(rows2)+"\n")
    rows4 = ["System,Bond,d(A),rho(e/A3),lap(e/A5),rho_BCP(e/bohr3),G(au),V(au),H(au),src"]
    cord = {}
    for r in sorted(IC, key=lambda r: -r[1]):
        dA, rb, lb, a, b, d1 = r[:6]
        g, v, h = gvh(rb, lb)
        ci = a if Z[a]==6 else b
        ii = a if Z[a]==53 else b
        cord.setdefault(ci, len(cord)+1)
        rows4.append(f"{label},C{cord[ci]}-I{Iorder[ii]},{dA:.3f},{rb/BOHR3:.3f},{lb/BOHR5:+.3f},{rb:.4f},{g:.4f},{v:.4f},{h:.4f},{r[7]}")
    open(f'{OUT}/tables/BCP_IC_full_442.csv','a').write("\n".join(rows4)+"\n")
    rows3 = ["System,Bond,d(A),rho(e/A3),lap(e/A5),rho_BCP(e/bohr3),G(au),V(au),H(au),src"]
    for r in sorted(MI, key=lambda r: -r[1]):
        dA, rb, lb, a, b, d1 = r[:6]
        g, v, h = gvh(rb, lb)
        y = a if Z[a]==39 else b
        ii = a if Z[a]==53 else b
        rows3.append(f"{label},Y-I{Iorder[ii]},{dA:.3f},{rb/BOHR3:.3f},{lb/BOHR5:+.3f},{rb:.4f},{g:.4f},{v:.4f},{h:.4f},{r[8]}")
    open(f'{OUT}/tables/BCP_MI_native_442.csv','a').write("\n".join(rows3)+"\n")
    print(f'  II={len(II)} (shell-filtered {nshell}) IO={len(IO)} YI={len(MI)} other={nother}')
    if agree:
        print(f'  line-min vs rho_BCP (bonded I-I): within 3%: {sum(1 for x in agree if x<3)}/{len(agree)}, within 10%: {sum(1 for x in agree if x<10)}/{len(agree)}, median {np.median(agree):.1f}%')

open(f'{OUT}/tables/BCP_II_full_442.csv','w').close()
open(f'{OUT}/tables/BCP_IO_full_442.csv','w').close()
open(f'{OUT}/tables/BCP_IC_full_442.csv','w').close()
open(f'{OUT}/tables/BCP_MI_native_442.csv','w').close()
for label, d, tag in SYS:
    analyse(label, d, tag)
print('DONE')
