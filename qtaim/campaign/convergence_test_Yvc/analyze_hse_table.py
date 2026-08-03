#!/usr/bin/env python3
"""HSE06 Table 1: activated binding energies per I2, chain(polyiodide) vs molecular.
DeltaE_bind = (E_loaded - E_empty - 4*E_gas)/4 * 1312.75  [kJ/mol per I2].
poly-vs-mol stabilization = (E_mol - E_poly)/4*1312.75 (reference-cancelling, the
clean HSE number). Y poly/mol are the pre-existing hse_poly/mol_yact; empty+gas+Gd/Tb
from the hse_bcd queue."""
import re,os
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
RYKJ=1312.75; N=4
def E(tag):
    f=f"{D}/{tag}.out"
    if not os.path.exists(f): return None
    m=re.findall(r"^!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f,errors="ignore").read(),re.M)
    return float(m[-1]) if m else None
gas=E("hse_i2_vac")
print("="*60+"\nHSE06 activated binding table (kJ/mol per I2)\n"+"="*60)
print(f"gas-I2 HSE ref = {gas} Ry" + ("   [MISSING]" if gas is None else ""))
print(f"\n{'metal':6s}{'chain/I2':>11s}{'mol/I2':>11s}{'poly stab':>11s}")
print("-"*42)
rows=[]
for m,tp,tm,te in [("Y","hse_poly_yact","hse_mol_yact","hse_empty_yact"),
                   ("Gd","hse_poly_gdact","hse_mol_gdact","hse_empty_gdact"),
                   ("Tb","hse_poly_tbact","hse_mol_tbact","hse_empty_tbact")]:
    ep,em,ee=E(tp),E(tm),E(te)
    miss=[k for k,v in [('gas',gas),('poly',ep),('mol',em),('empty',ee)] if v is None]
    if miss: print(f"{m:6s}  (missing: {','.join(miss)})"); continue
    dp=(ep-ee-N*gas)/N*RYKJ; dm=(em-ee-N*gas)/N*RYKJ; st=(em-ep)/N*RYKJ
    print(f"{m:6s}{dp:>11.1f}{dm:>11.1f}{st:>11.1f}")
    rows.append((m,ep,em,ee))
if rows:
    print("\nraw HSE energies (Ry):")
    for m,ep,em,ee in rows: print(f"  {m}: poly={ep:.5f} mol={em:.5f} empty={ee:.5f}")
