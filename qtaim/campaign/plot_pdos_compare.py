import re, numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Liberation Sans","font.size":10,"axes.labelsize":11,
    "xtick.labelsize":9,"ytick.labelsize":9,"legend.fontsize":8.5,"legend.frameon":False,
    "axes.linewidth":1.0,"pdf.fonttype":42})
B="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
def Ef(s):
    m=re.findall(r"the Fermi energy is\s*(-?\d+\.\d+)",open(f"{B}/{s}/scf_loaded.out",errors="ignore").read())
    return float(m[-1]) if m else 0.0
def parse_in(s):
    t=open(f"{B}/{s}/scf_loaded.in").read()
    cellA=np.array([list(map(float,re.findall(r"-?\d+\.\d+",l))) for l in t.split("CELL_PARAMETERS")[1].splitlines()[1:4]]); cellT=cellA.T
    apos=[]
    for l in t.split("ATOMIC_POSITIONS")[1].split("K_POINTS")[0].splitlines():
        m=re.match(r"\s*([A-Z][a-z]?)\s+(-?\d.\d+)\s+(-?\d.\d+)\s+(-?\d.\d+)",l)
        if m: apos.append((m.group(1),np.array([float(m.group(2)),float(m.group(3)),float(m.group(4))])))
    return cellT,apos
def coo_O(s):
    cellT,apos=parse_in(s); P={i:f for i,(_,f) in enumerate(apos)}
    def dist(i,fc): df=apos[i][1]-fc; df-=np.round(df); return np.linalg.norm(cellT@df)
    bt=open(f"{B}/{s}/bcp.out",errors="ignore").read(); oset=set()
    for l in bt.splitlines():
        m=re.match(r"\s*\d+\s+\w+\s+\(3,-1\)\s+\S+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+\d+\s+\w+\s+(-?\d+\.\d+E[+-]\d+)\s+\S+\s+(-?\d+\.\d+E[+-]\d+)",l)
        if not m: continue
        fc=np.array([float(m.group(1)),float(m.group(2)),float(m.group(3))]); rho=float(m.group(4)); lap=float(m.group(5))
        if rho<0.003 or rho>0.06 or abs(lap)>0.5: continue
        ds=sorted((dist(i,fc),i) for i in P)[:2]; sp=[apos[ds[0][1]][0],apos[ds[1][1]][0]]
        if sp.count("I")!=1 or "O" not in sp: continue
        a,b=ds[0][1],ds[1][1]; ax=cellT@((apos[b][1]-apos[a][1])-np.round(apos[b][1]-apos[a][1])); L=np.linalg.norm(ax)
        if L>3.5: continue
        v=cellT@((fc-apos[a][1])-np.round(fc-apos[a][1])); along=v@(ax/L); perp=np.linalg.norm(v-along*(ax/L))
        if not(0.1<along<0.9*L) or perp>0.5: continue
        oa=b if apos[a][0]=="I" else a
        if apos[oa][0]=="O": oset.add(oa)
    return sorted(oset)
def sumwfc(s,atoms,wfcid="#2(p)"):
    apos=parse_in(s)[1]; d=f"{B}/{s}"; E=None; tot=None
    for i in atoms:
        sp=apos[i][0]
        fn=[f for f in os.listdir(d) if re.match(rf"pdos\.pdos_atm#{i+1}\({sp}\)_wfc{re.escape(wfcid)}",f)]
        if not fn: continue
        a=np.loadtxt(f"{d}/{fn[0]}",comments="#")
        if E is None: E=a[:,0]; tot=np.zeros(len(a))
        tot+=a[:,1]
    return E,tot

data={}
for s in ["Y_str","Y_act"]:
    ef=Ef(s); apos=parse_in(s)[1]
    Iat=[i for i,(sp,_) in enumerate(apos) if sp=="I"]; Oat=coo_O(s)
    EI,I5p=sumwfc(s,Iat); EO,O2p=sumwfc(s,Oat)
    data[s]={"ef":ef,"EI":EI-ef,"I5p":I5p,"EO":EO-ef,"O2p_perO":O2p/len(Oat),"nO":len(Oat)}

# --- plot: two panels (left = I-5p, right = O-2p per-oxygen), both with strained vs activated ---
fig,axes=plt.subplots(1,2,figsize=(10,4))
# left: I-5p
ax=axes[0]
ax.plot(data["Y_str"]["EI"],data["Y_str"]["I5p"],color="#9467bd",lw=1.5,ls="-",label="I 5p — strained")
ax.plot(data["Y_act"]["EI"],data["Y_act"]["I5p"],color="#9467bd",lw=1.5,ls="--",label="I 5p — activated")
ax.axvline(0,color="k",lw=0.5,ls=":"); ax.fill_between(data["Y_str"]["EI"],0,data["Y_str"]["I5p"],
    where=data["Y_str"]["EI"]<0,color="#9467bd",alpha=0.15)
ax.set_xlim(-6,5); ax.set_xlabel("E − E$_F$, eV"); ax.set_ylabel("PDOS, states/eV"); ax.set_title("I 5p")
ax.legend(loc="upper right"); ax.grid(ls="--",alpha=0.3)
# right: O-2p per-oxygen
ax=axes[1]
ax.plot(data["Y_str"]["EO"],data["Y_str"]["O2p_perO"],color="#d62728",lw=1.5,ls="-",label=f"O 2p (COO) — strained, {data['Y_str']['nO']} O")
ax.plot(data["Y_act"]["EO"],data["Y_act"]["O2p_perO"],color="#d62728",lw=1.5,ls="--",label=f"O 2p (COO) — activated, {data['Y_act']['nO']} O")
ax.axvline(0,color="k",lw=0.5,ls=":"); ax.fill_between(data["Y_str"]["EO"],0,data["Y_str"]["O2p_perO"],
    where=data["Y_str"]["EO"]<0,color="#d62728",alpha=0.15)
ax.set_xlim(-6,5); ax.set_xlabel("E − E$_F$, eV"); ax.set_ylabel("PDOS per O, states/eV"); ax.set_title("O 2p (carboxylate, per-oxygen)")
ax.legend(loc="upper right"); ax.grid(ls="--",alpha=0.3)
fig.suptitle("Y-MOF: strained vs activated — PDOS comparison",fontsize=12,y=1.01)
fig.tight_layout()
fig.savefig(f"{B}/Y_pdos_compare.png",dpi=600,bbox_inches="tight")
fig.savefig(f"{B}/Y_pdos_compare.pdf",bbox_inches="tight")
print("wrote Y_pdos_compare.png/.pdf (I-5p left, O-2p per-O right; solid=strained, dashed=activated)")
