#!/usr/bin/env python3
"""Write the corrected QTAIM I-I BCP table as a .docx (Table Grid, thin black borders),
matching the style of the existing QTAIM_bond_table.docx. NEW file only."""
from docx import Document
from docx.shared import Pt
from docx.oxml.ns import qn
from docx.enum.table import WD_ALIGN_VERTICAL
import re,os
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
SYS=[("Y_str","Y-MOF strained"),("Tb_str","Tb-MOF strained"),("Gd_str","Gd-MOF strained"),
     ("Y_act","Y-MOF activated"),("Tb_act","Tb-MOF activated"),("Gd_act","Gd-MOF activated")]
BOHR=0.5291772109

def parse_corr(sys):
    ii=[];io=[];sec=None
    for l in open(f"{BASE}/{sys}/{sys}_bcp_corrected.txt"):
        if "CORRECTED I-I" in l: sec="ii"; continue
        if "CORRECTED I-O" in l: sec="io"; continue
        c=l.split()
        if l.startswith("#") or not c or not c[0].startswith("I"): continue
        if len(c)>=7:
            rec=(c[0],float(c[1]),float(c[2]),float(c[3]),float(c[4]),float(c[5]),c[6])
            (ii if sec=="ii" else io).append(rec)
    return ii,io

doc=Document()
st=doc.styles["Normal"]; st.font.name="Liberation Sans"; st.font.size=Pt(10)

def set_grid(table):
    tbl=table._tbl; tblPr=tbl.tblPr
    borders=tblPr.find(qn("w:tblBorders"))
    if borders is None:
        from docx.oxml import OxmlElement
        borders=OxmlElement("w:tblBorders"); tblPr.append(borders)
    for e in ("top","left","bottom","right","insideH","insideV"):
        from docx.oxml import OxmlElement
        el=borders.find(qn(f"w:{e}"))
        if el is None:
            el=OxmlElement(f"w:{e}"); borders.append(el)
        el.set(qn("w:val"),"single"); el.set(qn("w:sz"),"4")
        el.set(qn("w:space"),"0"); el.set(qn("w:color"),"000000")

p=doc.add_paragraph(); r=p.add_run("QTAIM I–I bond critical points — corrected (PAW/shell artifacts removed)")
r.bold=True; r.font.size=Pt(11)
doc.add_paragraph("ρ_BCP in e·Å⁻³, ∇²ρ in e·Å⁻⁵. Re-located on the clean valence field (plot_num=0) "
                  "with ρ17, pooled; CPs with same-atom gradient termini (intra-atomic 5p-shell), "
                  "near-nuclear residuals, or ρ outside the bonding window were excluded. ρ decreases "
                  "monotonically with d for every system; the bond-path profile agrees "
                  "(BCP_changelog.txt). EML local bond energies are omitted (the ½|V_BCP| relation is "
                  "calibrated on light-atom/H-bond systems and is unreliable for heavy iodine).")

for sys,name in SYS:
    ii,io=parse_corr(sys); ii=sorted(ii,key=lambda x:x[1])
    doc.add_paragraph(""); h=doc.add_paragraph(); hr=h.add_run(name); hr.bold=True
    t=doc.add_table(rows=1,cols=5); t.style="Table Grid"; set_grid(t)
    cells=t.rows[0].cells
    for c,txt in zip(cells,["I–I bond","d (Å)","ρ_BCP (e·Å⁻³)","∇²ρ (e·Å⁻⁵)","type"]):
        c.paragraphs[0].add_run(txt).bold=True
    for pair,d,rho,lap,rhoA,lapA,typ in ii:
        bt="covalent" if lap<-0.02 else ("closed-shell" if lap>0.02 else "borderline")
        row=t.add_row().cells
        vals=[pair,f"{d:.2f}",f"{rhoA:.3f}",f"{lapA:+.3f}",bt]
        for c,v in zip(row,vals):
            c.paragraphs[0].add_run(v)
            c.vertical_alignment=WD_ALIGN_VERTICAL.CENTER

doc.save(f"{BASE}/QTAIM_bond_table_corrected.docx")
print("wrote QTAIM_bond_table_corrected.docx")

# ---------------- I-O corrected table (same rigour as I-I) ----------------
doc2=Document()
s2=doc2.styles["Normal"]; s2.font.name="Liberation Sans"; s2.font.size=Pt(10)
p=doc2.add_paragraph(); r=p.add_run("QTAIM I···O bond critical points — corrected (same rigour as I–I)")
r.bold=True; r.font.size=Pt(11)
doc2.add_paragraph("ρ_BCP in e·Å⁻³, ∇²ρ in e·Å⁻⁵. Located on clean valence field (plot_num=0)+ρ17, pooled; "
   "intra-atomic 5p-shell CPs (ρ≈0.10, same-atom termini), near-nuclear residuals (∇²ρ≈+1.5), and "
   "out-of-window CPs excluded. Each BCP confirmed by an independent bond-path profile (agreement <1%). "
   "Contacts are predominantly at carboxylate (–COO⁻) oxygens; Gd-activated has one μ₃-OH contact.")
for sys,name in SYS:
    ii,io=parse_corr(sys); io=sorted(io,key=lambda x:x[1])
    if not io: continue
    doc2.add_paragraph(""); h=doc2.add_paragraph(); hr=h.add_run(name); hr.bold=True
    t=doc2.add_table(rows=1,cols=4); t.style="Table Grid"; set_grid(t)
    for c,txt in zip(t.rows[0].cells,["I···O bond","d (Å)","ρ_BCP (e·Å⁻³)","∇²ρ (e·Å⁻⁵)"]):
        c.paragraphs[0].add_run(txt).bold=True
    for pair,d,rho,lap,rhoA,lapA,typ in io:
        row=t.add_row().cells
        for c,v in zip(row,[pair,f"{d:.2f}",f"{rhoA:.3f}",f"{lapA:+.3f}"]):
            c.paragraphs[0].add_run(v); c.vertical_alignment=WD_ALIGN_VERTICAL.CENTER
doc2.save(f"{BASE}/IO_bond_table_corrected.docx")
print("wrote IO_bond_table_corrected.docx")
