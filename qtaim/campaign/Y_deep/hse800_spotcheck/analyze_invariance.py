#!/usr/bin/env python3
"""Final invariance analysis for I2@Y-MOFvc: PBE vs HSE06.
 (a) Precise pointwise valence-density diff (plot_num=0 cubes).
 (b) Per-bond I-I / I-O BCP rho comparison (plot_num=17 cubes), mapped via fractional MIC.
Both cubes are on the identical 150x144x320 grid, identical geometry."""
import re, numpy as np
SP="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/hse800_spotcheck"
YD="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep"
BOHR=0.5291772109
def parse_cube(path):
    L=open(path).readlines()
    nat=int(L[2].split()[0])
    nr=np.array([int(L[i].split()[0]) for i in (3,4,5)])
    step=np.array([[float(x) for x in L[i].split()[1:4]] for i in (3,4,5)])  # bohr/voxel, along each axis
    cellB=nr[:,None]*step                      # lattice vectors (bohr), rows = a1,a2,a3
    sym=[]; cartB=[]
    for i in range(6,6+nat):
        c=L[i].split(); sym.append(c[0] if not c[0].isdigit() else {1:"H",6:"C",8:"O",39:"Y",53:"I"}.get(int(c[0]),"X"))
        cartB.append([float(c[2]),float(c[3]),float(c[4])])
    cartB=np.array(cartB)                       # bohr
    frac=cartB @ np.linalg.inv(cellB)           # atom fractional coords
    vals=re.findall(r"[-+]?\d+\.\d+[eE][-+]?\d+|[-+]?\d+\.\d+", "".join(L[6+nat:]))
    data=np.array([float(v) for v in vals]).reshape(nr)
    return sym,frac,cellB,nr,data
def parse_bcp(outpath):
    rows=[]
    for l in open(outpath,errors="ignore"):
        m=re.match(r"\s*\d+\s+\S+\s+\(3,-1\)\s+bond\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+\d+\s+\S+\s+(-?\d+\.\d+E[+-]\d+)\s+\S+\s+(-?\d+\.\d+E[+-]\d+)",l)
        if m:
            rows.append((np.array([float(m.group(1)),float(m.group(2)),float(m.group(3))]),float(m.group(4)),float(m.group(5))))
    return rows
def bond_lbl(fc,frac,sym):
    df=frac-fc; df=df-np.round(df)        # fractional MIC
    d=np.linalg.norm(df,axis=1)           # fractional distance (cell-normalized, fine for ranking)
    i=np.argsort(d)[:2]
    return "".join(sorted(sym[i[0]]+sym[i[1]])), d[i[0]]+d[i[1]]
# ---------- (a) density diff ----------
print("="*70+"\n(a) VALENCE DENSITY DIFF: PBE val.cube vs HSE hse_val.cube (plot_num=0)\n"+"="*70)
_,_,_,nrP,pbe=parse_cube(f"{YD}/val.cube")
_,_,_,nrH,hse=parse_cube(f"{SP}/hse_val.cube")
assert tuple(nrP)==tuple(nrH)
d=pbe-hse
rmsP=np.sqrt(np.mean(pbe**2)); rmsD=np.sqrt(np.mean(d**2))
print(f"  grid {tuple(nrP)}, {pbe.size} voxels")
print(f"  RMS(rho)               = {rmsP:.6f} e/bohr3")
print(f"  RMS(|rho_PBE-rho_HSE|) = {rmsD:.3e} e/bohr3   ({100*rmsD/rmsP:.4f}% of RMS)")
print(f"  max|rho_PBE-rho_HSE|   = {np.abs(d).max():.3e} e/bohr3")
for lo,hi,lab in [(0.003,0.02,"I-O / halogen-bond region (0.003-0.02)"),
                  (0.03,0.09,"I-I polyiodide region (0.03-0.09)"),
                  (0.005,0.12,"all bonding BCP region (0.005-0.12)")]:
    m=(pbe>=lo)&(pbe<hi)
    if m.sum(): print(f"  {lab}: n={m.sum():6d}  mean|rPBE-rHSE|/rho = {np.mean(np.abs(d[m])/pbe[m])*100:.3f}%")
# ---------- (b) per-bond BCP comparison ----------
print("\n"+"="*70+"\n(b) BCP rho_BCP per bond: PBE rho17 vs HSE rho17 (plot_num=17)\n"+"="*70)
sym,frac,cellB,_,_=parse_cube(f"{YD}/rho17.cube")
P=parse_bcp(f"{SP}/bcp_PBE17.out"); H=parse_bcp(f"{SP}/bcp_HSE17.out")
def classified(rows):
    out={'II':[],'IO':[]}
    for fc,rho,lap in rows:
        lbl,s=bond_lbl(fc,frac,sym)
        if lbl=='II': out['II'].append((fc,rho,lap))
        elif lbl in('IO','OI'): out['IO'].append((fc,rho,lap))
    return out
pc=classified(P); hc=classified(H)
print(f"  PBE: {len(P)} BCPs total -> I-I={len(pc['II'])}, I-O={len(pc['IO'])}")
print(f"  HSE: {len(H)} BCPs total -> I-I={len(hc['II'])}, I-O={len(hc['IO'])}")
def match_and_compare(a,b,lab):
    # match each PBE BCP to nearest HSE BCP in fractional space
    a=sorted(a,key=lambda r:-r[1]); b=np.array([r[0] for r in b])
    rows=[]
    for fc,rho,lap in a:
        if len(b)==0: continue
        df=b-fc; df=df-np.round(df); j=np.argmin(np.linalg.norm(df,axis=1))
        # nearest HSE BCP rho
        hrho=[r[1] for r in (a*1)] # placeholder
        rows.append((rho,lap))
    return a
print("\n  --- I-I polyiodide BCPs (sorted by rho, PBE vs nearest HSE) ---")
pi=sorted(pc['II'],key=lambda r:-r[1]); hi=hc['II']
print(f"  {'rho_PBE':>9} {'rho_HSE':>9} {'dPBE-HSE':>9} {'%diff':>7}   {'lap_PBE':>9} {'lap_HSE':>9}")
for fc,rho,lap in pi:
    # nearest HSE I-I BCP
    if hi:
        arr=np.array([h[0] for h in hi]); df=arr-fc; df=df-np.round(df); j=np.argmin(np.linalg.norm(df,axis=1))
        hr=hi[j][1]; hl=hi[j][2]
        print(f"  {rho:9.4f} {hr:9.4f} {rho-hr:9.4f} {100*(rho-hr)/rho:7.2f}   {lap:9.2e} {hl:9.2e}")
print("\n  --- I-O halogen-bond BCPs (sorted by rho, PBE vs nearest HSE) ---")
po=sorted(pc['IO'],key=lambda r:-r[1]); ho=hc['IO']
for fc,rho,lap in po:
    if ho:
        arr=np.array([h[0] for h in ho]); df=arr-fc; df=df-np.round(df); j=np.argmin(np.linalg.norm(df,axis=1))
        hr=ho[j][1]; hl=ho[j][2]
        print(f"  {rho:9.4f} {hr:9.4f} {rho-hr:9.4f} {100*(rho-hr)/rho:7.2f}   {lap:9.2e} {hl:9.2e}")
# summary
if pi and hi:
    pr=np.array([r[1] for r in pi]);
    matched=[]
    for fc,rho,lap in pi:
        arr=np.array([h[0] for h in hi]); df=arr-fc; df=df-np.round(df); j=np.argmin(np.linalg.norm(df,axis=1)); matched.append(hi[j][1])
    matched=np.array(matched)
    print(f"\n  I-I mean rho: PBE={pr.mean():.4f}  HSE={matched.mean():.4f}  mean|%diff|={np.mean(np.abs(pr-matched)/pr)*100:.2f}%")
