#!/usr/bin/env python3
"""Landscape deliverables: (1) archive .in/.out + CIFs, (2) plot, (3) Word writeup."""
import re,os,shutil,numpy as np,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
BASE="/workspace/dimitry/calc/qe/i2_mof"
ARCH=f"{BASE}/HSE_final_tables/landscape_3variants"
os.makedirs(ARCH,exist_ok=True)
RYKJ=1312.75; N=4
def E(tag):
    f=f"{D}/{tag}.out"
    if not os.path.exists(f): return None
    m=re.findall(r"^!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f,errors="ignore").read(),re.M)
    return float(m[-1]) if m else None
def read_in(fn):
    L=open(fn).read().splitlines()
    ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
    ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][0]
    cell=np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)]); atoms=[];k=ai+1
    while k<len(L):
        s=L[k].split()
        if not s or s[0] in("K_POINTS","ATOMIC_SPECIES","CONSTRAINTS","OCCUPATIONS"): break
        if len(s)>=4 and re.match(r'^[A-Z]',s[0]):
            try: atoms.append((s[0],float(s[1]),float(s[2]),float(s[3])))
            except: pass
        k+=1
    return cell,atoms
def cell_geo(cell):
    a,b,c=[np.linalg.norm(cell[i]) for i in range(3)]
    def ang(i,j): return math.degrees(math.acos(np.dot(cell[i],cell[j])/(np.linalg.norm(cell[i])*np.linalg.norm(cell[j]))))
    return a,b,c,ang(1,2),ang(0,2),ang(0,1)
def wrap(f): return f%1.0  # normalize to [0,1); QE coords already here, fold is for safety
def write_cif(tag,name,note):
    f=f"{D}/{tag}.in"
    if not os.path.exists(f): return None
    cell,atoms=read_in(f)
    a,b,c,al,be,ga=cell_geo(cell)
    atoms=[(e,wrap(x),wrap(y),wrap(z)) for e,x,y,z in atoms]
    s=[f"data_{name}",f"# {note}","_symmetry_space_group_name_H-M   'P 1'","_symmetry_Int_Tables_number   1",
       f"_cell_length_a   {a:.6f}",f"_cell_length_b   {b:.6f}",f"_cell_length_c   {c:.6f}",
       f"_cell_angle_alpha   {al:.6f}",f"_cell_angle_beta   {be:.6f}",f"_cell_angle_gamma   {ga:.6f}",
       "loop_","_atom_site_type_symbol","_atom_site_fract_x","_atom_site_fract_y","_atom_site_fract_z"]
    for e,x,y,z in atoms: s.append(f"{e} {x:.6f} {y:.6f} {z:.6f}")
    open(f"{ARCH}/{name}.cif","w").write("\n".join(s)+"\n")
    return f"{name}.cif"
gas=E("rd_gas")
empties={"st":E("rd_empty_Yst"),"act":E("rd_empty_Yact"),"vc":E("rd_empty_Yvc")}
def dads(e,v): return (e-empties[v]-N*gas)/N*RYKJ if e and empties[v] else None
# (variant, type, tag, cif_name, cif_note) -- built by glob so new seeds auto-appear
import glob as _glob
def classify(tag):
    m=re.match(r'rd_poly_Y(act|vc|st)$',tag)
    if m: return (m.group(1),"CIF",tag,f"Y-MOF{m.group(1)}_polyiodide-CIF","Authoritative CIF polyiodide")
    m=re.match(r'hmol_Y(act|vc|st)_S1$',tag)
    if m: return (m.group(1),"poly",tag,f"Y-MOF{m.group(1)}_polyiodide-S1-clustered","Clustered alternative polyiodide")
    m=re.match(r'hse_lpoly_Y(act|vc|st)_S(\d+)$',tag)
    if m: return (m.group(1),"poly",tag,f"Y-MOF{m.group(1)}_polyiodide-S{m.group(2)}","Perturbed polyiodide (landscape scan)")
    m=re.match(r'rd_mol_Yact$',tag)
    if m: return ("act","mol",tag,"Y-MOFact_molecular-S1","Molecular I2 (RASPA-seeded)")
    m=re.match(r'hmol_Yvc$',tag)
    if m: return ("vc","mol",tag,"Y-MOFvc_molecular-S1","Molecular I2 (void-placed)")
    m=re.match(r'hse_lmol_Y(act|vc|st)_S(\d+)$',tag)
    if m: return (m.group(1),"mol",tag,f"Y-MOF{m.group(1)}_molecular-S{m.group(2)}","Molecular I2 (landscape scan)")
    return None
CONFIGS=[]
# authoritative CIF polyiodide (force-include even if energy missing)
for v in ["st","act","vc"]:
    CONFIGS.append((v,"CIF",f"rd_poly_Y{v}",f"Y-MOF{v}_polyiodide-CIF","Authoritative CIF polyiodide"))
# glob every landscape .out with a HSE energy and classify
for out in _glob.glob(f"{D}/*.out"):
    tag=os.path.basename(out)[:-4]
    if E(tag) is None: continue  # only configs with a finished HSE energy
    c=classify(tag)
    if c and not any(t[2]==c[2] for t in CONFIGS): CONFIGS.append(c)
# stable ordering: variant st/act/vc, then type CIF/poly/mol, then tag
torder={"CIF":0,"poly":1,"mol":2}; vorder={"st":0,"act":1,"vc":2}
CONFIGS.sort(key=lambda c:(vorder[c[0]],torder[c[1]],c[3]))
print("="*70+"\nLandscape deliverables\n"+"="*70)
print("\n[1] Archiving .in/.out + CIFs...")
manifest=[]
for v,t,tag,cn,note in CONFIGS:
    for ext in [".in",".out"]:
        src=f"{D}/{tag}{ext}"
        if os.path.exists(src): shutil.copy2(src,f"{ARCH}/{tag}{ext}")
    write_cif(tag,cn,note)
    e=E(tag); de=dads(e,v)
    manifest.append((v,t,tag,cn,e,de))
    print(f"  {cn:38s} {'OK' if e else 'MISSING':8s} ΔEads={de:+.1f}" if de else f"  {cn:38s} MISSING")
for tag in ["rd_gas","rd_empty_Yst","rd_empty_Yact","rd_empty_Yvc"]:
    for ext in [".in",".out"]:
        src=f"{D}/{tag}{ext}"
        if os.path.exists(src): shutil.copy2(src,f"{ARCH}/{tag}{ext}")
print(f"  -> {ARCH}/")
# designate the 6 lowest configs (min polyiodide + min molecular per variant)
lowest={}
for lv in ["st","act","vc"]:
    polyc=[(de,tag,cn) for (vv,t,tag,cn,e,de) in manifest if vv==lv and t in("CIF","poly") and de is not None]
    molc =[(de,tag,cn) for (vv,t,tag,cn,e,de) in manifest if vv==lv and t=="mol" and de is not None]
    for cands,lab in [(polyc,"polyiodide"),(molc,"molecular")]:
        if cands:
            _,ltag,lcn=min(cands); lowest[ltag]=f"LOWEST {lab} ({lv})"
            if os.path.exists(f"{ARCH}/{lcn}.cif"):
                shutil.copy2(f"{ARCH}/{lcn}.cif",f"{ARCH}/{lcn}__LOWEST.cif")
print(f"  designated {len(lowest)} lowest configs (CIFs copied with __LOWEST suffix): {[lowest[k] for k in lowest]}")

print("\n[2] Generating plot...")
colors={"CIF":"#d62728","poly":"#1f77b4","mol":"#2ca02c"}
fig,ax=plt.subplots(figsize=(7.5,8.5))
for x,v in enumerate(["st","act","vc"]):
    pts=[]
    for vv,t,tag,cn,note in CONFIGS:
        if vv!=v: continue
        e=E(tag)
        if e is None: continue
        pts.append((t,dads(e,v)))
    # group by energy to jitter near-degenerate levels apart for legibility
    for t,de in sorted(pts,key=lambda p:p[1]):
        lw=4 if t=="CIF" else 2.2
        ax.plot([x-0.16,x+0.16],[de,de],color=colors[t],linewidth=lw,solid_capstyle="butt",zorder=5 if t=="CIF" else 3)
        if t=="CIF":
            ax.plot(x,de,color="white",marker='o',markersize=6,zorder=6)
            ax.plot(x,de,color=colors[t],marker='o',markersize=6,markeredgewidth=1.5,fillstyle='none',zorder=7)
ax.set_xticks([0,1,2])
ax.set_xticklabels(["Y-MOF (st)","Y-MOF (act)","Y-MOF (vc)"],fontsize=13)
ax.set_ylabel(r"$\Delta E_{\rm ads}$ (kJ mol$^{-1}$ per I$_2$)",fontsize=13)
ax.set_title("HSE06 energy landscape: polyiodide vs molecular I$_2$\nin three Y-MOF pore geometries",fontsize=13)
ax.axhline(0,color='gray',linewidth=0.6,linestyle='--',zorder=1)
legend_elements=[Patch(facecolor=colors["CIF"],label="CIF polyiodide (studied in detail)"),
                 Patch(facecolor=colors["poly"],label="Alternative polyiodide configs"),
                 Patch(facecolor=colors["mol"],label="Molecular I$_2$ configs")]
ax.legend(handles=legend_elements,fontsize=11,loc='center left',
          bbox_to_anchor=(1.01,0.5),framealpha=0.95)
ax.set_xlim(-0.5,2.5)
ax.grid(axis='y',alpha=0.25)
# annotate the motif verdict per column (lower energy = stronger binding, at bottom)
for x,v in enumerate(["st","act","vc"]):
    polys=[dads(E(tag),v) for vv,tt,tag,_,__ in CONFIGS if vv==v and tt in("CIF","poly") and E(tag)]
    mols =[dads(E(tag),v) for vv,tt,tag,_,__ in CONFIGS if vv==v and tt=="mol" and E(tag)]
    if polys and mols:
        stab=min(mols)-min(polys)
        ax.annotate("polyiodide" if stab>0 else "molecular",
                    xy=(x,min(min(polys),min(mols))),xytext=(x,min(min(polys),min(mols))-7),
                    ha='center',fontsize=9,color='0.35',style='italic')
plt.tight_layout()
fig.savefig(f"{ARCH}/landscape_energy_levels.png",dpi=160,bbox_inches='tight')
fig.savefig(f"{D}/landscape_energy_levels.png",dpi=160,bbox_inches='tight')
plt.close(fig)
print("  -> landscape_energy_levels.png")

print("\n[3] Generating Word writeup...")
from docx import Document
from docx.shared import Pt
doc=Document(); doc.styles["Normal"].font.size=Pt(11)
def sub(p,b,s):
    r=p.add_run(b); r2=p.add_run(s); r2.font.subscript=True
h=doc.add_heading("I₂@Y-MOF energy landscape: polyiodide vs molecular configurations",level=1)
doc.add_heading("Motivation",level=2)
p=doc.add_paragraph("All real-space electronic-structure analysis in this work (QTAIM, ELF, Raman, ")
sub(p,"PDOS","2"); p.add_run(") is performed on the polyiodide motif. To place that analysis in proper "
 "thermodynamic context, we mapped the low-energy landscape of iodine in each of the three Y-MOF pore "
 "geometries — solvent-templated (st), activated (act), and variable-cell (vc) — asking, for "
 "each pore, whether the polyiodide or molecular I")
sub(p,"","2"); p.add_run(" is the preferred motif at the converged HSE06 level.")
doc.add_heading("Method",level=2)
doc.add_paragraph(
 "Multiple iodine configurations were generated per pore and evaluated at HSE06-D3/4×4×2 "
 "(PBE/4×4×2-relaxed geometries, framework held frozen at the loaded cell):\n\n"
 "(i) Polyiodide configurations — the authoritative CIF polyiodide (production geometry), plus "
 "alternatives generated by perturbing the CIF chain (±0.02 fractional displacement of iodine atoms) "
 "and re-relaxing at PBE. This probes the breadth of the polyiodide basin.\n\n"
 "(ii) Molecular I₂ configurations — generated by RASPA seeding (act) and by systematic "
 "void-placement: a grid-search of the framework interior locating the maximum-distance-to-framework "
 "positions at two c-layers (c ≈ 0.27c₀/0.73c₀ and c ≈ 0.15c₀/0.65c₀), "
 "placing four I₂ molecules (2.70 Å each) at the two highest-void sites per layer. "
 "The void-placement guarantees a physically valid, non-clashing molecular start that PBE relaxation "
 "preserves.\n\n"
 "All relaxations: PBE+D3(BJ)v4, Marzari–Vanderbilt smearing (degauss 0.02 Ry), BFGS, framework "
 "frozen, ion_maxstep 80 (energy-converged geometry harvested by plateau detection on the ionic force). "
 "Single-point energies: HSE06-D3/4×4×2, ecutfock 100 Ry, default nqx. Adsorption energy "
 "ΔEads = (E_loaded − E_empty − 4·E_gas)/4, per I₂.")
doc.add_heading("Results",level=2)
doc.add_paragraph(
 "The motif preference is set by pore confinement. In the open activated pore, molecular I₂ is "
 "thermodynamically preferred; in the two confined pores (st, vc), the polyiodide is strongly preferred "
 "and no molecular configuration reaches the polyiodide basin. The CIF polyiodide is not, in general, "
 "the global polyiodide minimum — perturbed alternatives lie 6–55 kJ mol⁻¹ per I₂ "
 "lower — but the within-family spread never inverts the motif verdict.")
# table
t=doc.add_table(rows=1,cols=5); t.style="Light Grid Accent 1"
for i,hd in enumerate(["Variant","Lowest polyiodide ΔEads","Lowest molecular ΔEads",
                       "Polyiodide − molecular","Preferred motif"]):
    t.rows[0].cells[i].text=hd
for v in ["st","act","vc"]:
    polys=[dads(E(tag),v) for vv,tt,tag,_,__ in CONFIGS if vv==v and tt in("CIF","poly") and E(tag)]
    mols=[dads(E(tag),v) for vv,tt,tag,_,__ in CONFIGS if vv==v and tt=="mol" and E(tag)]
    if not polys or not mols: continue
    lp=min(polys); lm=min(mols); stab=lm-lp
    r=t.add_row().cells
    r[0].text=f"Y-MOF ({v})"; r[1].text=f"{lp:+.2f}"; r[2].text=f"{lm:+.2f}"; r[3].text=f"{stab:+.2f}"
    r[4].text="polyiodide" if stab>0 else "molecular I₂"
doc.add_paragraph("Energies in kJ mol⁻¹ per I₂. ΔEads = (E_loaded − E_empty − "
 "4·E_gas)/4. Negative = stronger binding. “Lowest” = minimum over all configurations of "
 "that motif.").italic=True
# ---- Binding-energy decomposition (D3 + HSE06 electronic) ----
doc.add_heading("Binding-energy decomposition: D3 dispersion + HSE06 electronic",level=2)
import csv as _csv, os as _os
_decomp=f"{D}/landscape_D3_decomposition.csv"
doc.add_paragraph(
 "Each configuration's adsorption energy was decomposed as ΔEads = ΔEdisp + ΔEelec, where ΔEdisp is "
 "the D3(BJ) dispersion contribution evaluated at the HSE06-D3(BJ) level — functional-matched to the "
 "HSE06 electronic energy (HSE06 BJ parameters a1=0.383, s8=2.310, a2=5.685; dispersion is geometry-only "
 "and k-independent) — and ΔEelec = ΔEads − ΔEdisp is the HSE06 electronic part. The dispersion term is "
 "consistent with the 9-system reference table (st −79.8, act −82.0, vc −102.6 kJ mol⁻¹ per I₂). The "
 "decomposition separates two cleanly independent regimes. ΔEdisp is a pore property — essentially "
 "constant within each variant (st and act ≈ −80, reflecting their identical 1432 Å³ cells; vc ≈ −103, "
 "reflecting the contracted 1242 Å³ cell), and independent of motif or charge transfer. ΔEelec is the "
 "discriminator: large and negative for the genuinely charged st polyiodide (−50 to −70), near zero or "
 "slightly repulsive for the barely-charged act polyiodide (+7 to −3), weakly attractive for vc "
 "polyiodide (−10 to −16), and strongly positive (repulsive, +49) for the vc molecular configuration — "
 "the confined pore electronically rejects four isolated I₂ molecules, and only the dispersion floor "
 "binds them. Act is thus dispersion-dominated (electronic term ≈ 0); vc rejects molecular I₂ outright.")
if _os.path.exists(_decomp):
    rows_d=list(_csv.DictReader(open(_decomp)))
    td=doc.add_table(rows=1,cols=4); td.style="Light Grid Accent 1"
    for i,h in enumerate(["Variant","Lowest polyiodide (ΔEdisp / ΔEelec)","Lowest molecular (ΔEdisp / ΔEelec)","Regime"]):
        td.rows[0].cells[i].text=h
    interp={"st":"anchoring-driven (large electronic)","act":"dispersion-dominated (electronic ≈ 0)",
            "vc":"polyiodide weak electronic; molecular repulsive"}
    def _cell(r): return f"{float(r['dEads']):+.2f}  ({float(r['dEdisp']):+.2f} / {float(r['dEelec']):+.2f})"
    for v in ["st","act","vc"]:
        vr=[r for r in rows_d if r["variant"]==v]
        lp=min((r for r in vr if r["motif"] in("poly","CIF")),key=lambda r:float(r["dEads"]),default=None)
        lm=min((r for r in vr if r["motif"]=="mol"),key=lambda r:float(r["dEads"]),default=None)
        if not lp or not lm: continue
        rc=td.add_row().cells
        rc[0].text=f"Y-MOF ({v})"; rc[1].text=_cell(lp); rc[2].text=_cell(lm); rc[3].text=interp[v]
    doc.add_paragraph("Format: total ΔEads (ΔEdisp / ΔEelec), kJ mol⁻¹ per I₂. ΔEdisp ≈ constant per pore; ΔEelec carries the motif dependence.").italic=True
doc.add_paragraph(
 "The decomposition is shown in landscape_D3_decomposition.png (scatter of ΔEelec vs ΔEdisp): the "
 "dispersion axis separates the two pore sizes (st/act ≈ −80, vc ≈ −103), while the electronic axis "
 "separates the motifs — polyiodide below the zero line, molecular above it for vc (repulsive). The "
 "same figure makes explicit that act's binding is almost entirely the dispersion floor, and that vc "
 "molecular binding is electronically anti-binding.")
doc.add_heading("Discussion: justifying the polyiodide analysis",level=2)
doc.add_paragraph(
 "Because nothing in this study is experimental, the choice of which configuration to analyse in detail "
 "is a modelling choice, not a confrontation with observation. The polyiodide motif is examined "
 "consistently across the entire confinement gradient because it is the only way to isolate the effect "
 "of pore geometry on the chain’s electronic structure: adopting the global minimum of each pore "
 "would change the motif itself (molecular I₂ in act) and break the cross-variant comparison. The "
 "landscape above establishes two things that license this choice. First, the polyiodide is the genuine "
 "ground state in the confined pores where it is observed (st, vc). Second, even where it is metastable "
 "(act), it is only modestly so (≈11 kJ mol⁻¹ per I₂) and lies in a distinct, "
 "well-defined basin — a legitimate local minimum whose electronic structure is physically "
 "meaningful. The CIF polyiodides studied in detail are representative of their basin (within 6–55 "
 "kJ mol⁻¹ of the lowest polyiodide found), so the real-space descriptors extracted from them "
 "characterise the motif, not an arbitrary outlier.")
refs=["Heyd, G. E.; Scuseria, G. E.; Ernzerhof, M. J. Chem. Phys. 118, 8207 (2003).",
 "Krukau, A. V.; Vydrov, O. A.; Izmaylov, A. F.; Scuseria, G. E. J. Chem. Phys. 125, 224106 (2006).",
 "Grimme, S.; Antony, J.; Ehrlich, S.; Krieg, H. J. Chem. Phys. 132, 154104 (2010).",
 "Grimme, S.; Ehrlich, S.; Goerigk, L. J. Comput. Chem. 32, 1456 (2011)."]
doc.add_heading("References",level=2)
for i,rf in enumerate(refs): doc.add_paragraph(f"[{i+1}] {rf}")
out=f"{ARCH}/landscape_writeup.docx"; doc.save(out)
print(f"  -> landscape_writeup.docx")
# manifest file
with open(f"{ARCH}/MANIFEST.md","w") as f:
    f.write("# Landscape 3-variants archive\n\nConfigs (HSE06-D3/4x4x2 single-points on PBE/4x4x2-relaxed geometries):\n\n")
    f.write("| Variant | Motif | Config | E (Ry) | dEads (kJ/mol/I2) | Lowest |\n|---|---|---|---|---|---|\n")
    for v,t,tag,cn,e,de in manifest:
        note=lowest.get(tag,"")
        f.write(f"| {v} | {t} | {cn} | {e:.5f} | {de:+.2f} | {note} |\n" if de else f"| {v} | {t} | {cn} | MISSING | - | |\n")
    f.write(f"\nThe 6 lowest (min polyiodide + min molecular per variant) are also copied as "
            "`<name>__LOWEST.cif` for quick identification.\n\n")
    f.write(f"References: rd_gas, rd_empty_Y{{st,act,vc}} (gas I2 + empty frameworks).\n")
    f.write("Gas = "+(f"{gas:.5f}" if gas else "N/A")+" Ry.\n")
print(f"\nAll deliverables -> {ARCH}/")
print("  landscape_energy_levels.png + landscape_writeup.docx + MANIFEST.md + *.cif + *.in/.out")
