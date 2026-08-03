#!/usr/bin/env python3
"""Decompose landscape ΔEads into ΔEdisp(D3) + ΔEelec(HSE06) for all configs.
D3 from d3only_<tag>.out (PBE Gamma SP, geometry-only); HSE energy from <tag>.out."""
import re,os,glob,numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
RYKJ=1312.75; N=4
def E(tag):
    f=f"{D}/{tag}.out"
    m=re.findall(r"^!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f,errors="ignore").read(),re.M)
    return float(m[-1]) if m else None
def D3(tag):
    # D3 is geometry-only & k-independent: try the d3only Gamma SP first, then the
    # source PBE-relax .out (same final geometry; relax prints D3 every ionic step).
    hmmap={"hmol_Yvc":"mmol_Yvc","hmol_Yvc_S1":"mseed_Yvc_S1","hmol_Yst_S1":"mseed_Yst_S1"}
    srcs=[f"d3only_{tag}",
          (tag[4:] if tag.startswith("hse_") else None),
          hmmap.get(tag)]
    for src in srcs:
        if not src: continue
        f=f"{D}/{src}.out"
        if not os.path.exists(f): continue
        m=re.findall(r"DFT-D3 Dispersion\s*=\s*(-?\d+\.\d+)",open(f,errors="ignore").read())
        if m: return float(m[-1])
    return None
def classify(tag):
    m=re.match(r'rd_poly_Y(act|vc|st)$',tag);
    if m: return (m.group(1),"CIF")
    if re.match(r'hse_lpoly|hse_lmol',tag):
        v=re.search(r'Y(act|vc|st)',tag).group(1); return (v,"mol" if "lmol" in tag else "poly")
    if re.match(r'hmol_Y(act|vc|st)_S1$',tag): return (re.search(r'Y(act|vc|st)',tag).group(1),"poly")
    if re.match(r'rd_mol_Yact$',tag): return ("act","mol")
    if re.match(r'hmol_Yvc$',tag): return ("vc","mol")
    return None
gas=E("rd_gas"); d3gas=D3("rd_gas")
EMPT={"st":E("rd_empty_Yst"),"act":E("rd_empty_Yact"),"vc":E("rd_empty_Yvc")}
DEM={"st":D3("rd_empty_Yst"),"act":D3("rd_empty_Yact"),"vc":D3("rd_empty_Yvc")}
print(f"refs: gas E={gas} D3={d3gas}  | empty D3 st/act/vc = {DEM}")
rows=[]
for out in glob.glob(f"{D}/*.out"):
    tag=os.path.basename(out)[:-4]
    c=classify(tag)
    if c is None: continue
    e=E(tag); d3=D3(tag)
    if e is None or d3 is None: continue
    v,mt=c; dads=(e-EMPT[v]-N*gas)/N*RYKJ; disp=(d3-DEM[v]-N*d3gas)/N*RYKJ; elec=dads-disp
    rows.append({"sys":f"Y-MOF{v}","variant":v,"motif":mt,"tag":tag,"dEads":dads,"dEdisp":disp,"dEelec":elec})
# validate vs hardcoded 3-config values
HARDCODED={"st":-64.01,"act":-64.48,"vc":-83.46}
print("\n=== validate ΔEdisp(rd_poly) vs hardcoded 9-system table ===")
for v in ["st","act","vc"]:
    r=next((r for r in rows if r["tag"]==f"rd_poly_Y{v}"),None)
    if r: print(f"  Y-{v}: ΔEdisp={r['dEdisp']:.2f}  (hardcoded {HARDCODED[v]:.2f})  ΔEelec={r['dEelec']:.2f}")
vorder={"st":0,"act":1,"vc":2}; mord={"CIF":0,"poly":1,"mol":2}
rows.sort(key=lambda r:(vorder[r["variant"]],mord[r["motif"]],r["dEads"]))
# table
with open(f"{D}/landscape_D3_decomposition.csv","w") as f:
    f.write("variant,motif,tag,dEads,dEdisp,dEelec\n")
    for r in rows: f.write(f"{r['variant']},{r['motif']},{r['tag']},{r['dEads']:.2f},{r['dEdisp']:.2f},{r['dEelec']:.2f}\n")
print(f"\nwrote landscape_D3_decomposition.csv ({len(rows)} configs)")
# plot: scatter ΔEelec (y) vs ΔEdisp (x) — both terms visible, motif=color, variant=marker
from matplotlib.lines import Line2D as L2
fig,ax=plt.subplots(figsize=(9,8))
mcol={"poly":"#1f77b4","CIF":"#1f77b4","mol":"#d62728"}  # polyiodide (incl CIF)=blue, molecular=red
vmark={"st":"o","act":"s","vc":"^"}
xx=np.array([-92,-58])
for c in (-50,-100,-130):
    ax.plot(xx,c-xx,color='0.75',ls='--',lw=0.9,zorder=1)
    ax.text(-58.5,c+58.5,f"$\\Sigma$={c}",fontsize=8,color='0.55',va='center')
for r in rows:
    cif=r["motif"]=="CIF"
    ax.scatter(r["dEdisp"],r["dEelec"],c=mcol[r["motif"]],marker=('*' if cif else vmark[r["variant"]]),
               s=(220 if cif else 70),edgecolor=('goldenrod' if cif else 'k'),linewidth=(2.2 if cif else 0.5),
               zorder=5 if cif else 3)
ax.axhline(0,color='k',lw=1)
ax.axvline(-73.5,color='0.65',ls='-',lw=0.8)
ax.text(-64.5,-73.5,"st/act pore",fontsize=8,color='0.5',ha='center',va='top')
ax.text(-83.5,-73.5,"vc pore",fontsize=8,color='0.5',ha='center',va='top')
ax.set_xlabel(r"$\Delta E_{\rm disp}$  (D3 dispersion, kJ mol$^{-1}$ per I$_2$)",fontsize=12)
ax.set_ylabel(r"$\Delta E_{\rm elec}$  (HSE06 electronic, kJ mol$^{-1}$ per I$_2$)",fontsize=12)
ax.set_title("Landscape decomposition: dispersion (pore-size) vs electronic (motif)\n"
             "blue=polyiodide, red=molecular; ○ st  □ act  △ vc  ★ = CIF (studied)",fontsize=11)
ax.text(-83,55,"vc molecular\n(repulsive electronic)",color='#d62728',fontsize=9,ha='center')
ax.text(-63.5,-9,"act polyiodide\n(dispersion-dominated)",color='#1f77b4',fontsize=9,ha='center')
ax.text(-63.5,-72,"st polyiodide\n(strong electronic)",color='#1f77b4',fontsize=9,ha='center')
ax.grid(alpha=0.2,zorder=0)
ax.set_xlim(-92,-58)
handles=[L2([0],[0],marker='o',color='w',markerfacecolor='#1f77b4',label='polyiodide',markersize=10),
         L2([0],[0],marker='o',color='w',markerfacecolor='#d62728',label='molecular',markersize=10),
         L2([0],[0],marker='*',color='w',markerfacecolor='#1f77b4',markeredgecolor='goldenrod',markeredgewidth=1.5,label='CIF (studied)',markersize=16),
         L2([0],[0],marker='o',color='w',markerfacecolor='0.6',label='st',markersize=9),
         L2([0],[0],marker='s',color='w',markerfacecolor='0.6',label='act',markersize=9),
         L2([0],[0],marker='^',color='w',markerfacecolor='0.6',label='vc',markersize=9)]
ax.legend(handles=handles,fontsize=9,ncol=2,loc='upper right',framealpha=0.95)
plt.tight_layout()
fig.savefig(f"{D}/landscape_D3_decomposition.png",dpi=150,bbox_inches='tight')
fig.savefig(f"/workspace/dimitry/calc/qe/i2_mof/HSE_final_tables/landscape_3variants/landscape_D3_decomposition.png",dpi=150,bbox_inches='tight')
print("wrote landscape_D3_decomposition.png")
# summary
print("\n=== per-variant lowest-poly vs lowest-mol decomposition ===")
for v in ["st","act","vc"]:
    for mt in ["poly","mol","CIF"]:
        rr=[r for r in rows if r["variant"]==v and (r["motif"]==mt or (mt=="poly" and r["motif"]=="CIF"))]
        if rr:
            lo=min(rr,key=lambda r:r["dEads"])
            print(f"  Y-{v} {mt:4s} lowest: ΔEads={lo['dEads']:.1f} = disp {lo['dEdisp']:.1f} + elec {lo['dEelec']:.1f}")
