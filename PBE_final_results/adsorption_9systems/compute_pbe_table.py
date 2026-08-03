#!/usr/bin/env python3
"""Generate the PBE-D3 9-system adsorption table. The SCF energies in the 'HSE' outputs
are PBE electronic (EXX was never added — the QE 7.5 PAW-EXX bug). So PBE-D3 ΔEads =
PBE_electronic + PBE-D3 dispersion = (HSE_reported - HSE06-D3) + PBE-D3.
All D3 values from the dftd3 driver (d3_pbe_vs_hse06.json)."""
import json
D3=json.load(open("/workspace/dimitry/calc/qe/i2_mof/HSE_final_tables/adsorption_9systems/d3_pbe_vs_hse06.json"))
RYKJ=1312.75; N=4
# HSE-reported ΔEads (PBE electronic + HSE06-D3) — from the binding table
HSE_DADS={("Y","st"):-130.08,("Y","act"):-69.64,("Y","vc"):-111.98,
          ("Gd","st"):-151.22,("Gd","act"):-89.12,("Gd","vc"):-117.18,
          ("Tb","st"):-131.74,("Tb","act"):-77.66,("Tb","vc"):-113.71}
def dD3(m,var,func):
    L=D3[f"rd_poly_{m}{var}"][func]; E=D3[f"rd_empty_{m}{var}"][func]; G=D3["rd_gas"][func]
    return (L-E-N*G)/N*RYKJ
print(f"{'System':16s} {'ΔEads_HSE':>10s} {'ΔEdisp_PBE':>11s} {'ΔEdisp_HSE06':>13s} {'ΔEads_PBE-D3':>13s} {'ΔEelec':>8s}")
print("-"*80)
results={}
for var in ["st","act","vc"]:
  for m in ["Gd","Tb","Y"]:
    hse=HSE_DADS[(m,var)]
    dp=dD3(m,var,"pbe"); dh=dD3(m,var,"hse06")
    pbe=hse+(dp-dh)  # swap D3 from HSE06 to PBE
    elec=pbe-dp      # electronic part (should be same as HSE elec)
    results[(m,var)]=pbe
    print(f"{m+'-MOF'+var:16s} {hse:10.2f} {dp:11.2f} {dh:13.2f} {pbe:13.2f} {elec:8.2f}")
  print()
# write CSV
import csv
with open("/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/adsorption_9systems/pbe_d3_table.csv","w") as f:
    w=csv.writer(f); w.writerow(["System","ΔEads_PBE-D3","ΔEdisp_PBE-D3","ΔEelec"])
    for var in ["st","act","vc"]:
        for m in ["Gd","Tb","Y"]:
            dp=dD3(m,var,"pbe"); pbe=results[(m,var)]
            w.writerow([f"{m}-MOF{var}",f"{pbe:.2f}",f"{dp:.2f}",f"{pbe-dp:.2f}"])
print("wrote pbe_d3_table.csv")
