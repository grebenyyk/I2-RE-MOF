#!/usr/bin/env python3
"""Write the corrected Table 1 (ΔEads = ΔED3 + ΔEHSE06, HSE06-D3) to CSV + DOCX,
for pasting into the manuscript. ΔEads unchanged; ΔED3 now HSE06-D3(BJ)."""
import json, os
from docx import Document
from docx.shared import Pt
A9 = "/workspace/dimitry/calc/qe/i2_mof/HSE_final_tables/adsorption_9systems"
d = json.load(open(f"{A9}/d3_pbe_vs_hse06.json"))
RYKJ = 1312.75; N = 4
DADS = {("Y","st"):-130.08,("Y","act"):-69.64,("Y","vc"):-111.98,
        ("Gd","st"):-151.22,("Gd","act"):-89.12,("Gd","vc"):-117.18,
        ("Tb","st"):-131.74,("Tb","act"):-77.66,("Tb","vc"):-113.71}
ORD=[("st",[("Gd","Gd-MOFst"),("Tb","Tb-MOFst"),("Y","Y-MOFst")]),
     ("act",[("Gd","Gd-MOFact"),("Tb","Tb-MOFact"),("Y","Y-MOFact")]),
     ("vc",[("Gd","Gd-MOFvc"),("Tb","Tb-MOFvc"),("Y","Y-MOFvc")])]
def dEd3(m,var):
    L=d[f"rd_poly_{m}{var}"]["hse06"]; E=d[f"rd_empty_{m}{var}"]["hse06"]; G=d["rd_gas"]["hse06"]
    return (L-E-N*G)/N*RYKJ
# CSV
with open(f"{A9}/table1_HSE06D3.csv","w") as f:
    f.write("System,ΔEads/I2,ΔED3/I2,ΔEHSE06/I2\n")
    for var,metals in ORD:
        for m,name in metals:
            da=DADS[(m,var)]; dd=dEd3(m,var); f.write(f"{name},{da:.2f},{dd:.2f},{da-dd:.2f}\n")
        f.write("\n")
# DOCX
doc=Document()
doc.add_heading("Table 1 (corrected, HSE06-D3)",level=2)
doc.add_paragraph("DFT-calculated electronic energies of adsorption ΔEads of iodine in RE-MOFs for the "
 "three constrained structural models. Each unit cell loaded with four I₂ (GCMC). ΔED3 is the Grimme "
 "D3(BJ) dispersion contribution functional-matched to HSE06 (HSE06 BJ params a1=0.383, s8=2.310, "
 "a2=5.685): ΔEads = ΔED3 + ΔEHSE06. All values kJ mol⁻¹ per I₂.")
doc.add_paragraph(
 "Footnote for the table: ΔED3 is evaluated with the HSE06-D3(BJ) parameters (s8=2.310, a1=0.383, "
 "a2=5.685; Grimme 2011, J. Comput. Chem. 32, 1456), functional-matched to the HSE06 electronic energy. "
 "These parameters are more attractive than the PBE-D3(BJ) set (s8=0.7875) by ~17 kJ mol⁻¹ per I₂, so "
 "the dispersion scale exceeds the PBE-D3 values commonly reported for iodine/MOF systems; the total "
 "ΔEads is independent of this decomposition.").italic=True
lab={"st":"Solvent-templated","act":"Activated","vc":"Variable-cell"}
for var,metals in ORD:
    doc.add_paragraph(lab[var]).bold=True
    t=doc.add_table(rows=1,cols=4); t.style="Light Grid Accent 1"
    for i,h in enumerate(["System","ΔEads/I₂","ΔED3/I₂","ΔEHSE06/I₂"]): t.rows[0].cells[i].text=h
    for m,name in metals:
        da=DADS[(m,var)]; dd=dEd3(m,var); rc=t.add_row().cells
        rc[0].text=name; rc[1].text=f"{da:.2f}"; rc[2].text=f"{dd:.2f}"; rc[3].text=f"{da-dd:.2f}"
    doc.add_paragraph("")
doc.save(f"{A9}/table1_HSE06D3.docx")
print("wrote table1_HSE06D3.csv + table1_HSE06D3.docx")
print(f"{'System':14s} {'dEads':>8s} {'dED3':>8s} {'dEHSE06':>8s}")
for var,metals in ORD:
    for m,name in metals:
        da=DADS[(m,var)]; dd=dEd3(m,var)
        print(f"{name:14s} {da:8.2f} {dd:8.2f} {da-dd:8.2f}")
    print()
