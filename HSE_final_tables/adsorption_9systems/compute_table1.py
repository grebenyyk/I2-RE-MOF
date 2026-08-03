#!/usr/bin/env python3
"""Recompute the 9-system Table 1 (ΔEads = ΔED3 + ΔEHSE06) with HSE06-D3 (correct)
vs PBE-D3 (current table). ΔEads is unchanged (HSE total energies already include
HSE06-D3); only the dispersion/electronic SPLIT changes."""
import json
A9 = "/workspace/dimitry/calc/qe/i2_mof/HSE_final_tables/adsorption_9systems"
d = json.load(open(f"{A9}/d3_pbe_vs_hse06.json"))
RYKJ = 1312.75; N = 4
# ΔEads (kJ/I2) from the locked HSE binding table (unchanged)
DADS = {("Y","st"):-130.08,("Y","act"):-69.64,("Y","vc"):-111.98,
        ("Gd","st"):-151.22,("Gd","act"):-89.12,("Gd","vc"):-117.18,
        ("Tb","st"):-131.74,("Tb","act"):-77.66,("Tb","vc"):-113.71}
def dEd3(metal,var,func):
    L = d[f"rd_poly_{metal}{var}"][func]; E = d[f"rd_empty_{metal}{var}"][func]; G = d["rd_gas"][func]
    return (L - E - N*G)/N * RYKJ
print(f"{'system':16s} {'dEads':>8s} | {'dED3_PBE':>9s} {'dEHSE_PBE':>10s} | {'dED3_HSE06':>11s} {'dEHSE06_HSE06':>13s}")
print("-"*78)
for var in ["st","act","vc"]:
    for metal in ["Y","Gd","Tb"]:
        da = DADS[(metal,var)]
        dp = dEd3(metal,var,"pbe"); dh = dEd3(metal,var,"hse06")
        print(f"{metal+'-MOF'+var:16s} {da:8.2f} | {dp:9.2f} {da-dp:10.2f} | {dh:11.2f} {da-dh:13.2f}")
    print()
