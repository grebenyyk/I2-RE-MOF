#!/usr/bin/env python3
"""Correct PBE-D3 landscape: E_PBE-D3 = E_HSE_reported - D3_HSE06 + D3_PBE.
This swaps the HSE06-D3 (baked into the HSE SCF totals) for PBE-D3,
giving genuine PBE-D3 total energies. Then ΔEads/ΔEdisp/ΔEelec are all
self-consistently PBE-D3."""
import re,os,glob,json,numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D as L2
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
OUT="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/landscape_3variants"
RYKJ=1312.75; N=4
HSE06=json.load(open(f"{D}/hse06d3_energies.json"))  # D3_HSE06 per config
PBE=json.load(open(f"{D}/pbe_d3_energies.json"))       # D3_PBE per config

def E(tag):
    f=f"{D}/{tag}.out"
    m=re.findall(r"^!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f,errors="ignore").read(),re.M)
    return float(m[-1]) if m else None
def d3h(tag):  # HSE06-D3
    v=HSE06.get(tag); return v if isinstance(v,(int,float)) else None
def d3p(tag):  # PBE-D3
    v=PBE.get(tag); return v if isinstance(v,(int,float)) else None
def Epbe(tag):  # E_PBE-D3 = E_HSE - D3_HSE06 + D3_PBE
    e=E(tag); h=d3h(tag); p=d3p(tag)
    if e is None or h is None or p is None: return None
    return e - h + p
def classify(tag):
    m=re.match(r'rd_poly_Y(act|vc|st)$',tag)
    if m: return (m.group(1),"CIF")
    if re.match(r'hse_lpoly|hse_lmol',tag):
        v=re.search(r'Y(act|vc|st)',tag).group(1); return (v,"mol" if "lmol" in tag else "poly")
    if re.match(r'hmol_Y(act|vc|st)_S1$',tag): return (re.search(r'Y(act|vc|st)',tag).group(1),"poly")
    if re.match(r'rd_mol_Yact$',tag): return ("act","mol")
    if re.match(r'hmol_Yvc$',tag): return ("vc","mol")
    return None

gas=Epbe("rd_gas"); d3gas_p=d3p("rd_gas")
EMPT={v:Epbe(f"rd_empty_Y{v}") for v in ["st","act","vc"]}
DEM={v:d3p(f"rd_empty_Y{v}") for v in ["st","act","vc"]}
print(f"refs: gas PBE-D3 E={gas} D3={d3gas_p}")
print(f"empty PBE-D3: {EMPT}")

rows=[]
for out in glob.glob(f"{D}/*.out"):
    tag=os.path.basename(out)[:-4]; c=classify(tag)
    if c is None: continue
    ep=Epbe(tag); d3=d3p(tag)
    if ep is None or d3 is None: continue
    v,mt=c
    dads=(ep-EMPT[v]-N*gas)/N*RYKJ; disp=(d3-DEM[v]-N*d3gas_p)/N*RYKJ; elec=dads-disp
    rows.append({"sys":f"Y-MOF{v}","variant":v,"motif":mt,"tag":tag,"dEads":dads,"dEdisp":disp,"dEelec":elec})

print(f"\n=== validate rd_poly_Yvc (expect ΔEads≈-92.89, ΔEdisp≈-83.46, ΔEelec≈-9.43) ===")
r=next((r for r in rows if r["tag"]=="rd_poly_Yvc"),None)
if r: print(f"  ΔEads={r['dEads']:.2f}  ΔEdisp={r['dEdisp']:.2f}  ΔEelec={r['dEelec']:.2f}")

vorder={"st":0,"act":1,"vc":2}; mord={"CIF":0,"poly":1,"mol":2}
rows.sort(key=lambda r:(vorder[r["variant"]],mord[r["motif"]],r["dEads"]))
with open(f"{OUT}/landscape_D3_decomposition.csv","w") as f:
    f.write("variant,motif,tag,dEads,dEdisp,dEelec\n")
    for r in rows: f.write(f"{r['variant']},{r['motif']},{r['tag']},{r['dEads']:.2f},{r['dEdisp']:.2f},{r['dEelec']:.2f}\n")
print(f"\nwrote landscape_D3_decomposition.csv ({len(rows)} configs) [correct PBE-D3]")

# plot
fig,ax=plt.subplots(figsize=(9,8))
mcol={"poly":"#1f77b4","CIF":"#1f77b4","mol":"#d62728"}; vmark={"st":"o","act":"s","vc":"^"}
xx=np.array([-92,-58])
for c in (-30,-70,-100):
    ax.plot(xx,c-xx,color='0.75',ls='--',lw=0.9,zorder=1)
    ax.text(-58.5,c+58.5,f"$\\Sigma$={c}",fontsize=8,color='0.55',va='center')
for r in rows:
    cif=r["motif"]=="CIF"
    ax.scatter(r["dEdisp"],r["dEelec"],c=mcol[r["motif"]],marker=('*' if cif else vmark[r["variant"]]),
               s=(220 if cif else 70),edgecolor=('goldenrod' if cif else 'k'),linewidth=(2.2 if cif else 0.5),zorder=5 if cif else 3)
ax.axhline(0,color='k',lw=1)
ax.set_xlabel(r"$\Delta E_{\rm disp}$  (PBE-D3, kJ mol$^{-1}$ per I$_2$)",fontsize=12)
ax.set_ylabel(r"$\Delta E_{\rm elec}$  (PBE electronic, kJ mol$^{-1}$ per I$_2$)",fontsize=12)
ax.set_title("Landscape decomposition (PBE-D3): dispersion vs electronic\n"
             "blue=polyiodide, red=molecular",fontsize=11)
ax.grid(alpha=0.2,zorder=0); ax.set_xlim(-92,-58)
handles=[L2([0],[0],marker='o',color='w',markerfacecolor='#1f77b4',label='polyiodide',markersize=10),
         L2([0],[0],marker='o',color='w',markerfacecolor='#d62728',label='molecular',markersize=10),
         L2([0],[0],marker='*',color='w',markerfacecolor='#1f77b4',markeredgecolor='goldenrod',markeredgewidth=1.5,label='CIF',markersize=16),
         L2([0],[0],marker='o',color='w',markerfacecolor='0.6',label='st',markersize=9),
         L2([0],[0],marker='s',color='w',markerfacecolor='0.6',label='act',markersize=9),
         L2([0],[0],marker='^',color='w',markerfacecolor='0.6',label='vc',markersize=9)]
ax.legend(handles=handles,fontsize=9,ncol=2,loc='upper right',framealpha=0.95)
plt.tight_layout()
fig.savefig(f"{OUT}/landscape_D3_decomposition.png",dpi=150,bbox_inches='tight')
print("wrote landscape_D3_decomposition.png [correct PBE-D3]")

print("\n=== per-variant lowest-poly vs lowest-mol (PBE-D3) ===")
for v in ["st","act","vc"]:
    for mt in ["poly","mol","CIF"]:
        rr=[r for r in rows if r["variant"]==v and (r["motif"]==mt or (mt=="poly" and r["motif"]=="CIF"))]
        if rr:
            lo=min(rr,key=lambda r:r["dEads"])
            print(f"  Y-{v} {mt:4s} lowest: ΔEads={lo['dEads']:.1f} = disp {lo['dEdisp']:.1f} + elec {lo['dEelec']:.1f}")
