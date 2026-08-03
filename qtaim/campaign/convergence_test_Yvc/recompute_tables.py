#!/usr/bin/env python3
"""Recompute both binding tables from the 4x4x2 redo (self-consistent 50/400).
loaded/isubl @ 4x4x2 (new), empty @ Gamma (new), gasI2 + vc-relaxed empty @ Gamma (have).
D3 (geometry-only, k-independent) carried over from production -> DeltaE_D3 column unchanged."""
import re,os,numpy as np
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
RY=1312.75
def E(tag):
    f=os.path.join(D,tag+".out")
    if not os.path.exists(f): return None
    m=re.findall(r"^!.*total energy\s*=\s*(-?\d+\.\d+)",open(f).read(),re.M)
    return float(m[-1]) if m else None
GAS=float(re.findall(r"^!.*total energy\s*=\s*(-?\d+\.\d+)",open("/workspace/dimitry/calc/qe/i2_mof/i2_vac_ref/scf.out").read(),re.M)[-1])
# vc-relaxed empty (from vc_empty run) for MOF-relax
def Evcrel(v):
    f=f"/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/{v}/vc_empty.out"
    return float(re.findall(r"^!.*total energy\s*=\s*(-?\d+\.\d+)",open(f).read(),re.M)[-1])
DED3={"Y_str":-64.01,"Y_actc":-64.48,"Y_deep":-83.46,"Gd_str":-62.90,"Gd_act":-62.28,"Gd_deep":-84.50,"Tb_str":-63.12,"Tb_act":-63.51,"Tb_deep":-81.38}
VAR=[("Gd","Gd_str"),("Tb","Tb_str"),("Y","Y_str"),("Gd","Gd_act"),("Tb","Tb_act"),("Y","Y_actc"),("Gd","Gd_deep"),("Tb","Tb_deep"),("Y","Y_deep")]
GEO={"Gd_str":"st","Tb_str":"st","Y_str":"st","Gd_act":"act","Tb_act":"act","Y_actc":"act","Gd_deep":"vc","Tb_deep":"vc","Y_deep":"vc"}
def kJ(x): return x*RY
print("="*70+"\nTABLE 1  (4x4x2 / Gamma, 50/400 PBE+D3)  dE_ads = dE_D3 + dE_PBE\n"+"="*70)
print(f"{'system':14}{'dEads':>9}{'dED3':>9}{'dEPBE':>9}  (kJ/mol I2)")
for label,v in VAR:
    L=E(f"redo442_loaded_{v}"); F=E(f"redo442_empty_{v}")
    if L is None or F is None: print(f"{label+'-'+GEO[v]:14}  (missing)"); continue
    dads=kJ((L-F-4*GAS)/4); d3=DED3[v]; dpbe=dads-d3
    print(f"{'I2@'+label+'-MOF_'+GEO[v]:14}{dads:>9.2f}{d3:>9.2f}{dpbe:>9.2f}")
print("\n"+"="*70+"\nTABLE 2  chain-form / anchor / MOF-relax / dEint  (kJ/mol I2)\n"+"="*70)
print(f"{'system':14}{'chain':>9}{'anchor':>9}{'MOFrelax':>9}{'dEint':>9}")
T2=[("Y","Y_str"),("Y","Y_actc"),("Y","Y_deep"),("Tb","Tb_deep"),("Gd","Gd_deep")]
for label,v in T2:
    L=E(f"redo442_loaded_{v}"); I=E(f"redo442_isubl_{v}"); F=E(f"redo442_empty_{v}")
    if None in (L,I,F): print(f"{label+'-'+GEO[v]:14}  (missing)"); continue
    chain=kJ((I-4*GAS)/4); anchor=kJ((L-F-I)/4)
    if GEO[v]=="vc":
        relax=kJ((F-Evcrel(v))/4); dint=chain+anchor+relax; rstr=f"{relax:+.2f}"
    else:
        dint=chain+anchor; rstr="  -"
    print(f"{'I2@'+label+'-MOF_'+GEO[v]:14}{chain:>9.2f}{anchor:>9.2f}{rstr:>9}{dint:>9.2f}")
