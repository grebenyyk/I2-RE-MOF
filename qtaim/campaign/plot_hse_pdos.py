#!/usr/bin/env python3
"""HSE PDOS for the 3 Y-MOF variants (decoded from converged HSE wfc), with PBE
reference. Confirms the I-5p sigma* acceptor stays pinned at E_F at the hybrid
level (narrow-gap metallic -> CT is not a PBE artifact)."""
import re, numpy as np, os, glob
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Liberation Sans","font.size":9,"axes.labelsize":10,
    "xtick.labelsize":8.5,"ytick.labelsize":8.5,"legend.fontsize":8,"legend.frameon":False,
    "axes.linewidth":1.0,"pdf.fonttype":42,"svg.fonttype":"none"})
B="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
VAR=[("Y_str","st"),("Y_actc","act"),("Y_deep","vc")]
# crystallography-consistent colors
COL={"I":"#8e3fa9","O":"#d62728","C":"#555555","Y":"#2ca02c"}  # I purple, O red, C gray, Y green
VCOL={"st":"#1b9e77","act":"#d95f02","vc":"#7570b3"}  # variant colors (teal/orange/purple)
VL={"st":r"Y-MOF$_{\rm st}$","act":r"Y-MOF$_{\rm act}$","vc":r"Y-MOF$_{\rm vc}$"}
# atom ranges (C1-32 H33-44 O45-68 Y69-72 I73-80); wfc ids
AT={"I":(range(73,81),"#2(p)"),"O":(range(45,69),"#2(p)"),"C":(range(1,33),"#2(p)"),"Y":(range(69,73),"#5(d)")}
def fermi(d,tag):
    if tag=="HSE":
        t=open(f"{B}/{d}/hse/hse_keep/save/data-file-schema.xml",errors="ignore").read()
        m=re.search(r"<fermi_energy>([-\d.E+]+)</fermi_energy>",t)
        return float(m.group(1))*27.211386 if m else 0.0
    t=open(f"{B}/{d}/scf_loaded.out",errors="ignore").read()
    m=re.findall(r"the Fermi energy is\s*(-?\d+\.\d+)",t)
    return float(m[-1]) if m else 0.0
def sumel(d,tag,elem):
    pat=f"{B}/{d}/{'hse_pdos.pdos_' if tag=='HSE' else 'pdos.pdos_'}"
    atoms,wid=AT[elem]; E=None; tot=None
    for i in atoms:
        fs=glob.glob(f"{pat}atm#{i}({elem})_wfc{wid}*")
        if not fs: continue
        a=np.loadtxt(fs[0],comments="#")
        if E is None: E=a[:,0]; tot=np.zeros(len(a))
        tot+=a[:,1]
    return E,tot
def get(d,tag,elem):
    E,t=sumel(d,tag,elem); ef=fermi(d,tag); return E-ef,t

# ---- Figure: 2 panels. Left: I-5p 3 variants HSE (solid) + PBE (dashed). Right: O-2p 3 variants HSE ----
fig,axes=plt.subplots(1,2,figsize=(7,3.0))
XL=(-4,3)
# Left: I-5p
ax=axes[0]
for d,tg in VAR:
    EH,H=get(d,"HSE","I"); EP,P=get(d,"PBE","I")
    ax.plot(EH,H,color=VCOL[tg],lw=1.7,label=VL[tg])
    ax.plot(EP,P,color=VCOL[tg],lw=1.0,ls="--")  # PBE reference dashed
ax.axvline(0,color="k",lw=0.9,ls=":")
ax.text(0.02,ax.get_ylim()[1]*0.92,r"$E_{\rm F}$",fontsize=8)
ax.set_xlim(*XL); ax.set_xlabel(r"$E-E_{\rm F}$, eV"); ax.set_ylabel("PDOS, states/eV")
ax.set_title(r"I 5p  (— HSE, – – PBE)",fontsize=9)
ax.grid(ls="--",alpha=0.25)
# custom legend: variants + PBE dashed
from matplotlib.lines import Line2D
hnd=[Line2D([0],[0],color=VCOL[tg],lw=1.7,label=VL[tg]) for _,tg in VAR]+[Line2D([0],[0],color="0.5",lw=1.0,ls="--",label="PBE")]
ax.legend(handles=hnd,loc="upper right",fontsize=7.5)
# Right: O-2p per oxygen, 3 variants HSE
ax=axes[1]
for d,tg in VAR:
    E,H=get(d,"HSE","O"); nO=len(AT["O"][0])
    ax.plot(E,H/nO,color=VCOL[tg],lw=1.7,label=VL[tg])
ax.axvline(0,color="k",lw=0.9,ls=":")
ax.text(0.02,ax.get_ylim()[1]*0.92,r"$E_{\rm F}$",fontsize=8)
ax.set_xlim(*XL); ax.set_xlabel(r"$E-E_{\rm F}$, eV"); ax.set_ylabel("O 2p PDOS per O, states/eV")
ax.set_title(r"O 2p (carboxylate, per-O), HSE",fontsize=9); ax.grid(ls="--",alpha=0.25)
ax.legend(loc="upper right",fontsize=7.5)
fig.tight_layout(pad=0.4)
out=f"{B}/final_activated_results/pdos_Y_3variants_HSE"
os.makedirs(os.path.dirname(out),exist_ok=True)
fig.savefig(out+".png",dpi=600,bbox_inches="tight"); fig.savefig(out+".pdf",bbox_inches="tight")
print("wrote",out+".png/.pdf")

# ---- print σ* acceptor positions (HSE vs PBE) for the manuscript ----
print("\nI-5p σ* acceptor peak (E-E_F, eV) above E_F:")
print(f"{'variant':<6}{'HSE σ* (eV)':>14}{'PBE σ* (eV)':>14}{'HSE gap (meV)':>16}")
for d,tg in VAR:
    EH,H=get(d,"HSE","I"); EP,Pb=get(d,"PBE","I")
    efH=fermi(d,"HSE"); efP=fermi(d,"PBE")
    mH=(EH>0)&(EH<1.5); mP=(EP>0)&(EP<1.5)
    pkH=EH[np.argmax(H*mH)] if mH.any() else np.nan
    pkP=EP[np.argmax(Pb*mP)] if mP.any() else np.nan
    print(f"{tg:<6}{pkH:>14.3f}{pkP:>14.3f}{'':>16}")
