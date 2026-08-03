#!/usr/bin/env python3
"""Density-invariance spot-check (PBE vs HSE06) for I2@Y-MOFvc.
 (a) DIRECT pointwise density comparison on matched-grid plot_num=0 cubes (val.cube):
     RMS and max |rho_PBE - rho_HSE| relative to the density scale -> the headline.
 (b) critic2 BCP spectra: sorted rho_BCP values for I-I / I-O bonds, PBE vs HSE."""
import re, os, subprocess, numpy as np
SP="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/hse800_spotcheck"
YD="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep"
CRITIC="/workspace/dimitry/calc/qe/i2_mof/tools/critic2_run.sh"
def parse_cube(path):
    with open(path) as f: L=f.readlines()
    nat=int(L[2].split()[0])
    nr=[int(L[i].split()[0]) for i in (3,4,5)]
    cell=np.array([[float(x) for x in L[i].split()[1:4]] for i in (3,4,5)])*np.array(nr)[:,None]*0.529177 # step*nr->A
    atoms=[]
    for i in range(6,6+nat):
        c=L[i].split(); atoms.append((c[0],float(c[2])*0.529177,float(c[3])*0.529177,float(c[4])*0.529177))
    vals=re.findall(r"[-+]?\d+\.\d+[eE][-+]?\d+|[-+]?\d+\.\d+", "".join(L[6+nat:]))
    data=np.array([float(v) for v in vals]).reshape(nr)  # (n1,n2,n3), e/bohr3
    return cell,atoms,nr,data
# ---------- (a) direct density diff ----------
def direct_diff():
    pbe=parse_cube(f"{YD}/val.cube")
    hse=parse_cube(f"{SP}/hse_val.cube")
    assert pbe[2]==hse[2], f"grid mismatch {pbe[2]} vs {hse[2]}"
    a,b=pbe[3],hse[3]; d=a-b
    scale=np.sqrt(np.mean(a*a))   # RMS density (reference scale)
    print("=== (a) DIRECT pointwise density comparison: PBE val.cube vs HSE hse_val.cube ===")
    print(f"  grid {pbe[2]}, {a.size} voxels, units e/bohr3")
    print(f"  RMS(rho_PBE)            = {scale:.6f}")
    print(f"  RMS(rho_PBE - rho_HSE)  = {np.sqrt(np.mean(d*d)):.6f}   ({100*np.sqrt(np.mean(d*d))/scale:.2f}% of RMS)")
    print(f"  max |rho_PBE - rho_HSE| = {np.abs(d).max():.6f}")
    # where the biggest differences sit
    idx=np.unravel_index(np.abs(d).argmax(),d.shape)
    print(f"  max-diff voxel at {idx} (fractional {np.array(idx)/np.array(pbe[2])})")
    # BCP-region check: density in 0.005-0.12 e/bohr3 (bonding region) mean rel diff
    mask=(a>0.005)&(a<0.12)
    print(f"  bonding-region voxels (0.005<rho<0.12): mean|rPBE-rHSE|/rho = {np.mean(np.abs(d[mask])/a[mask])*100:.2f}%")
# ---------- (b) critic2 BCP spectra ----------
def bcp_spectrum(cube,tag):
    cri=f"crystal {cube}\nload {cube} id rho\nauto\n"
    open(f"{SP}/bcp_{tag}.cri","w").write(cri)
    subprocess.run(f"cd {SP} && {CRITIC} bcp_{tag}.cri > bcp_{tag}.out 2>&1",shell=True)
    t=open(f"{SP}/bcp_{tag}.out",errors="ignore").read()
    rh=[]
    for l in t.splitlines():
        m=re.match(r"\s*\d+\s+\w+\s+\(3,-1\)\s+\S+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+\d+\s+\w+\s+(-?\d+\.\d+E[+-]\d+)\s+\S+\s+(-?\d+\.\d+E[+-]\d+)",l)
        if m: rh.append(float(m.group(4)))
    return np.array(rh)
def bcp_compare():
    print("\n=== (b) critic2 BCP spectra: rho_BCP sorted desc (e/bohr3) ===")
    res={}
    for tag,path in [("PBE17",f"{YD}/rho17.cube"),("HSE17",f"{SP}/hse_rho17.cube"),
                     ("PBE0",f"{YD}/val.cube"),("HSE0",f"{SP}/hse_val.cube")]:
        if not os.path.exists(path): print(f"  MISSING {path}"); continue
        rh=bcp_spectrum(path,tag); res[tag]=rh
        s=np.sort(rh)[::-1]
        print(f"  {tag}: {len(rh)} BCPs | top15: "+" ".join(f"{v:.4f}" for v in s[:15]))
    # I-I vs I-O split: I-I BCPs are high-rho (>0.04), I-O are low (<0.015)
    print("\n  per-region mean rho_BCP (PBE vs HSE):")
    for lvl in ["17","0"]:
        for lo,hi,lab in [(0.04,9,"I-I (>0.04)"),(0.003,0.02,"I-O (0.003-0.02)")]:
            p=res.get("PBE"+lvl); h=res.get("HSE"+lvl)
            if p is None or h is None: continue
            pm=p[(p>=lo)&(p<hi)]; hm=h[(h>=lo)&(h<hi)]
            if len(pm): print(f"    plot{lvl} {lab}: PBE mean={pm.mean():.4f} (n={len(pm)})  HSE mean={hm.mean() if len(hm) else float('nan'):.4f} (n={len(hm)})")
if __name__=="__main__":
    direct_diff(); bcp_compare()
