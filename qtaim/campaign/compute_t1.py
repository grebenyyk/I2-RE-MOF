#!/usr/bin/env python3
"""Compute Table S2 (vc-relax binding corrections) at Table-1 level (50/400, dg0.01).
Fixed-cell leg = exact Table 1 values (campaign chain loaded + activate_empty).
VC leg = vc-relax_t1 (loaded) + scf_empty_vct2 (empty at vc cell). Fully consistent.
"""
import re,sys
RY2KJ=1312.75
CAMP="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
I2=-757.79989474; I2D=-0.00354048
IRT=4*I2; IRD=4*I2D
# Table 1 fixed-cell leg (verified to reproduce Table 1 exactly)
FIXED={
 "Gd":(-6865.15037671,-3833.54878663),  # activated_chain/loaded.out , activate_empty/relax.out
 "Tb":(-6925.19968642,-3893.61675089),  # activated_i2/loaded_ext.out , activate_empty/relax.out
 "Y" :(-5767.37949755,-2735.81928050),  # activated_i2_kickB/loaded.out, activate_empty/relax.out
}
DIR={"Gd":"Gd_act","Tb":"Tb_act","Y":"Y_actc"}

def grab(path):
    txt=open(path,errors="ignore").read()
    E=[float(x) for x in re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)",txt)]
    D=[float(x) for x in re.findall(r"DFT-D3 Dispersion\s*=\s*(-?\d+\.\d+)",txt)]
    C=re.findall(r"CELL_PARAMETERS[^\n]*\n((?:\s*-?[\d.eE+-]+[^\n]*\n){3})",txt)
    cells=[[ [float(z) for z in ln.split()[:3]] for ln in c.strip().split("\n")] for c in C]
    return E,D,cells
def vol(c):
    a,b,d=c
    return abs(a[0]*(b[1]*d[2]-b[2]*d[1])-a[1]*(b[0]*d[2]-b[2]*d[0])+a[2]*(b[0]*d[1]-b[1]*d[0]))

print(f"{'m':<3}{'E_fix':>9}{'E_vc':>10}{'strong':>8}{'E_D3(vc)':>9}{'E_elec(vc)':>11}{'  %dV':>7}")
out=[]
for m,d in DIR.items():
    lf,ef=FIXED[m]
    try:
        le,ld,lc=grab(f"{CAMP}/{d}/vc-relax_t1.out")
        eef,edf,_=grab(f"{CAMP}/{d}/scf_empty_vct2.out")   # empty at vc cell
    except Exception as ex:
        print(f"{m}: not ready ({ex})"); continue
    if not eef:
        print(f"{m}: empty_vct2 not converged yet"); continue
    lvc,evc=le[-1],eef[-1]
    bfix=(lf-ef-IRT)*RY2KJ/4
    bvc=(lvc-evc-IRT)*RY2KJ/4
    bvcd=(ld[-1]-edf[-1]-IRD)*RY2KJ/4
    vf=vol(lc[0]); vvc=vol(lc[-1]); pct=(vvc-vf)/vf*100
    out.append((m,bfix,bvc,bvc-bfix,bvcd,bvc-bvcd,pct,vf,vvc))
    print(f"{m:<3}{bfix:>9.2f}{bvc:>10.2f}{bvc-bfix:>8.2f}{bvcd:>9.2f}{bvc-bvcd:>11.2f}{pct:>6.1f}%")
if len(out)==3:
    import statistics
    print(f"\nmean strengthening = {statistics.mean([o[3] for o in out]):.2f} kJ/I2")
    print(f"check fixed == Table1: "+" ".join(f"{m}={o[1]:.2f}" for m,*_ in [('Gd',)]+[()]) if False else "")
    for m,bf,*_ in out:
        print(f"  {m} fixed {bf:.2f} (Table1: { {'Gd':-131.94,'Tb':-125.81,'Y':-118.36}[m] })")
