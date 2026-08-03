#!/usr/bin/env python3
"""HSE06 full 9-system binding table: st/act/vc x Y/Gd/Tb.
dEads/I2 = (E_loaded - E_empty - 4*E_gas)/4 * 1312.75  [kJ/mol per I2].
act already done; vc/st fill in as the two GPU queues complete."""
import re,os
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
RYKJ=1312.75; N=4
def E(tag):
    f=f"{D}/{tag}.out"
    if not os.path.exists(f): return None
    m=re.findall(r"^!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f,errors="ignore").read(),re.M)
    return float(m[-1]) if m else None
gas=E("hse_i2_vac")
METALS={"Y":"y","Gd":"gd","Tb":"tb"}
def tag(m,kind,var):  # kind=poly|empty, var=st|act|vc
    act=f"{m}act" if m!="y" else "yact"     # Y act uses 'yact', others '{m}act'
    suf={"act":act,"vc":f"{m}vc","st":f"{m}st"}[var]
    return f"hse_{kind}_{suf}"
print("="*60+"\nHSE06 full binding table (kJ/mol per I2)\n"+"="*60)
print(f"gas-I2 HSE = {gas} Ry"+("" if gas else "  [MISSING]"))
print(f"\n{'metal':5}{'st':>9}{'act':>9}{'vc':>9}")
print("-"*35)
done_all=True
for name,m in METALS.items():
    cells=[]
    for v in ["st","act","vc"]:
        L=E(tag(m,"poly",v)); F=E(tag(m,"empty",v))
        if None in (L,F,gas): cells.append(f"{'--':>9}"); done_all=False
        else: cells.append(f"{(L-F-N*gas)/N*RYKJ:>9.1f}")
    print(f"{name:5}{''.join(cells)}")
print(f"\nall 9 systems present: {done_all}")
