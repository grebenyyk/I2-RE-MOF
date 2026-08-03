#!/usr/bin/env python3
"""HSE06 REDO table (uniform settings). binding=(loaded-empty-4*gas)/4*1312.75.
loaded=rd_poly_{Metal}{st,act,vc} (4x4x2), empty=rd_empty_{...} (Gamma), gas=rd_gas.
dED3 = geometry-only (DED3, same as PBE table since CIF=production geom; HSE-D3≈PBE-D3<0.1kJ).
dEelec = dEads - dED3. All at degauss0.02/ef100/50-400/CIF geom (fixes Gd-act+vc-empty bugs,
mixed-mesh, Gamma-only-Fock, and the -13kJ gas ref)."""
import re,os
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
RYKJ=1312.75; N=4
DED3={("Gd","st"):-62.90,("Tb","st"):-63.12,("Y","st"):-64.01,
      ("Gd","act"):-62.28,("Tb","act"):-63.51,("Y","act"):-64.48,
      ("Gd","vc"):-84.50,("Tb","vc"):-81.38,("Y","vc"):-83.46}
def E(tag):
    f=f"{D}/{tag}.out"
    if not os.path.exists(f): return None
    m=re.findall(r"^!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f,errors="ignore").read(),re.M)
    return float(m[-1]) if m else None
gas=E("rd_gas")
print("="*64+"\nHSE06 REDO (uniform: 4x4x2 loaded, Gamma empty/gas, dg0.02, ef100, 50/400, CIF)\n"+"="*64)
print(f"gas (rd_gas) = {gas} Ry"+("" if gas else "  [MISSING]\n"))
print(f"\n{'System':16}{'dEads/I2':>10}{'dED3/I2':>10}{'dEelec/I2':>11}")
print("-"*48)
allok=gas is not None
MM=["Gd","Tb","Y"]
for v,vn in [("st","Solvent-templated"),("act","Activated"),("vc","Variable-cell")]:
    print(vn)
    for m in MM:
        L=E(f"rd_poly_{m}{v}"); F=E(f"rd_empty_{m}{v}")
        if None in (L,F,gas): print(f"  I2@{m}-MOF{v:<5}(missing: L={L is not None} F={F is not None})"); allok=False; continue
        dads=(L-F-N*gas)/N*RYKJ; d3=DED3[(m,v)]
        print(f"  I2@{m}-MOF{v:<5}{dads:>10.2f}{d3:>10.2f}{dads-d3:>11.2f}")
print(f"\nall 9 present: {allok}")
