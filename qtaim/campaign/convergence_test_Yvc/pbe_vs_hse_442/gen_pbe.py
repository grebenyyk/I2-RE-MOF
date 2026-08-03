#!/usr/bin/env python3
"""Generate PBE-D3 single-point inputs matching the HSE binding-table protocol,
by stripping the hybrid/EXX cards from the HSE inputs (same geometry/cutoff/k/D3)."""
import re,sys,os
D="/workspace/dimitry/calc/qe/i2_mof/HSE_final_tables/adsorption_9systems"
OUT="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/pbe_vs_hse_442"
os.makedirs(OUT,exist_ok=True)
STRIP=re.compile(r"^\s*(input_dft|exxdiv_treatment|x_gamma_extrapolation|exx_fraction|screening_parameter|ecutfock)\b",re.I)
for src,tag in [("rd_poly_Yvc.in","pbe_poly_Yvc"),("rd_empty_Yvc.in","pbe_empty_Yvc"),("rd_gas.in","pbe_gas")]:
    t=open(f"{D}/{src}").read().splitlines()
    out=[]
    for l in t:
        if STRIP.match(l): continue
        l=re.sub(r"prefix *= *'[^']*'",f"prefix = '{tag}'",l)
        out.append(l)
    open(f"{OUT}/{tag}.in","w").write("\n".join(out)+"\n")
    print("wrote",tag, [l.strip() for l in out if 'calculation' in l.lower() or 'K_POINTS' in l][:2])
