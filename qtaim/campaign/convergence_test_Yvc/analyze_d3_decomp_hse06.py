#!/usr/bin/env python3
"""Recompute landscape ΔEads = ΔEdisp + ΔEelec with HSE06-D3 (correct, functional-matched)
instead of PBE-D3. ΔEdisp now uses HSE06-D3(BJ) energies from hse06d3_energies.json
(evaluated via QE's bundled dftd3, func=hse06 version=4). Overwrites
landscape_D3_decomposition.csv + .png."""
import re,os,glob,json,numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D as L2
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
RYKJ=1312.75; N=4
HSE06=json.load(open(f"{D}/hse06d3_energies.json"))
def E(tag):
    f=f"{D}/{tag}.out"
    m=re.findall(r"^!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f,errors="ignore").read(),re.M)
    return float(m[-1]) if m else None
def D3(tag):
    v=HSE06.get(tag)
    return v if isinstance(v,(int,float)) else None
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
print(f"refs: gas E={gas} D3_HSE06={d3gas}  | empty D3 st/act/vc = {DEM}")
rows=[]
for out in glob.glob(f"{D}/*.out"):
    tag=os.path.basename(out)[:-4]; c=classify(tag)
    if c is None: continue
    e=E(tag); d3=D3(tag)
    if e is None or d3 is None: continue
    v,mt=c; dads=(e-EMPT[v]-N*gas)/N*RYKJ; disp=(d3-DEM[v]-N*d3gas)/N*RYKJ; elec=dads-disp
    rows.append({"sys":f"Y-MOF{v}","variant":v,"motif":mt,"tag":tag,"dEads":dads,"dEdisp":disp,"dEelec":elec})
# HSE06-D3 reference (rd_poly CIF) values for sanity print
print("\n=== ΔEdisp(rd_poly) HSE06-D3 per variant ===")
for v in ["st","act","vc"]:
    r=next((r for r in rows if r["tag"]==f"rd_poly_Y{v}"),None)
    if r: print(f"  Y-{v}: ΔEdisp_HSE06={r['dEdisp']:.2f}  ΔEelec={r['dEelec']:.2f}  ΔEads={r['dEads']:.2f}")
vorder={"st":0,"act":1,"vc":2}; mord={"CIF":0,"poly":1,"mol":2}
rows.sort(key=lambda r:(vorder[r["variant"]],mord[r["motif"]],r["dEads"]))
with open(f"{D}/landscape_D3_decomposition.csv","w") as f:
    f.write("variant,motif,tag,dEads,dEdisp,dEelec\n")
    for r in rows: f.write(f"{r['variant']},{r['motif']},{r['tag']},{r['dEads']:.2f},{r['dEdisp']:.2f},{r['dEelec']:.2f}\n")
print(f"\nwrote landscape_D3_decomposition.csv ({len(rows)} configs) [HSE06-D3]")
# plot
fig,ax=plt.subplots(figsize=(9,8))
mcol={"poly":"#1f77b4","CIF":"#1f77b4","mol":"#d62728"}; vmark={"st":"o","act":"s","vc":"^"}
xx=np.array([-112,-72])
for c in (-60,-100,-140):
    ax.plot(xx,c-xx,color='0.75',ls='--',lw=0.9,zorder=1)
    ax.text(-72,c+72,f"$\\Sigma$={c}",fontsize=8,color='0.55',va='center')
for r in rows:
    cif=r["motif"]=="CIF"
    ax.scatter(r["dEdisp"],r["dEelec"],c=mcol[r["motif"]],marker=('*' if cif else vmark[r["variant"]]),
               s=(220 if cif else 70),edgecolor=('goldenrod' if cif else 'k'),linewidth=(2.2 if cif else 0.5),
               zorder=5 if cif else 3)
ax.axhline(0,color='k',lw=1)
ax.axvline(-91,color='0.65',ls='-',lw=0.8)   # midpoint of st/act(~-80) and vc(~-103)
ax.text(-79,-92,"st/act pore",fontsize=8,color='0.5',ha='center',va='top')
ax.text(-103,-92,"vc pore",fontsize=8,color='0.5',ha='center',va='top')
ax.set_xlabel(r"$\Delta E_{\rm disp}$  (HSE06-D3, kJ mol$^{-1}$ per I$_2$)",fontsize=12)
ax.set_ylabel(r"$\Delta E_{\rm elec}$  (HSE06 electronic, kJ mol$^{-1}$ per I$_2$)",fontsize=12)
ax.set_title("Landscape decomposition (HSE06-D3): dispersion (pore-size) vs electronic (motif)\n"
             "blue=polyiodide, red=molecular; st  act  vc  = CIF",fontsize=11)
ax.grid(alpha=0.2,zorder=0)
ax.set_xlim(-112,-72)
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
print("wrote landscape_D3_decomposition.png [HSE06-D3]")
print("\n=== per-variant lowest-poly vs lowest-mol decomposition (HSE06-D3) ===")
for v in ["st","act","vc"]:
    for mt in ["poly","mol","CIF"]:
        rr=[r for r in rows if r["variant"]==v and (r["motif"]==mt or (mt=="poly" and r["motif"]=="CIF"))]
        if rr:
            lo=min(rr,key=lambda r:r["dEads"])
            print(f"  Y-{v} {mt:4s} lowest: ΔEads={lo['dEads']:.1f} = disp {lo['dEdisp']:.1f} + elec {lo['dEelec']:.1f}")
