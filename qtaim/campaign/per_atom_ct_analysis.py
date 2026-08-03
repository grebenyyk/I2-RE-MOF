#!/usr/bin/env python3
"""Per-atom charge transfer on the I-O bond network.
Merges BCP_IO_full.csv (I-O bond critical points) with critic2 ct.out (Bader-basin
deltarho integral per atom) -> charge on each bonded I and O. Outputs table + plot
of bond strength (rho_BCP) vs atom charge."""
import re,csv,os,numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
CAMP="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
BCP=f"{CAMP}/bcp_plots/multiwfn_out/BCP_IO_full.csv"
# system -> ct.out dir (Y_act uses the corrected Y_actc)
SYS2DIR={"Y-MOFst":"Y_str","Tb-MOFst":"Tb_str","Gd-MOFst":"Gd_str",
         "Y-MOFact":"Y_actc","Tb-MOFact":"Tb_act","Gd-MOFact":"Gd_act",
         "Y-MOFvc":"Y_deep","Tb-MOFvc":"Tb_deep","Gd-MOFvc":"Gd_deep"}
METAL={"Y":"Y","Tb":"Tb","Gd":"Gd"}; VAR={"st":"st","act":"act","vc":"vc"}
def parse_ct(sysname):
    """Return {label(e.g. 'I5'): dr} and per-element {elem: [dr,...]} from ct.out.
    $dr = last numeric column (4th integrable). Verified: I sum matches ct.txt."""
    d=f"{CAMP}/{SYS2DIR[sysname]}/ct.out"
    rows=[]
    pat=re.compile(r'^\s*(\d+)\s+\d+\s+\d+\s+([A-Z][a-z]?)\s+\d+\s+--\s+([0-9.+\-E]+)\s+([0-9.+\-E]+)\s+([0-9.+\-E]+)\s+([0-9.+\-E]+)')
    for l in open(d,errors="ignore"):
        m=pat.match(l)
        if m:
            iat=int(m.group(1)); el=m.group(2); dr=float(m.group(6))  # last col = $dr
            rows.append((iat,el,dr))
    rows.sort(key=lambda r:r[0])
    lbl={}; perel={}
    cnt={}
    for iat,el,dr in rows:
        cnt[el]=cnt.get(el,0)+1; label=f"{el}{cnt[el]}"; lbl[label]=dr
        perel.setdefault(el,[]).append(dr)
    return lbl,perel
# ---- load bonds ----
bonds=[]
with open(BCP) as f:
    for r in csv.DictReader(f):
        sys=r["System"]; bond=r["Bond"]
        m=re.match(r'(I\d+)-(O\d+)',bond)
        if not m: continue
        bonds.append({"sys":sys,"IL":m.group(1),"OL":m.group(2),
                      "d":float(r["d(A)"]),"rho":float(r["rho_BCP(e/bohr3)"])})
# ---- merge with per-atom charges ----
CT={}; AVG={}
for sys in SYS2DIR:
    lbl,perel=parse_ct(sys); CT[sys]=lbl
    AVG[sys]={el:(np.mean(v) if v else 0) for el,v in perel.items()}
    # sanity: I sum vs known ct.txt
    isum=sum(perel.get("I",[]))
    print(f"  {sys} ({SYS2DIR[sys]}): {len(lbl)} atoms, I sum={isum:+.3f}, avg I={AVG[sys].get('I',0):+.3f}, avg O={AVG[sys].get('O',0):+.3f}")
print()
# ---- build extended table ----
rows=[]
for b in bonds:
    Idr=CT[b["sys"]].get(b["IL"]); Odr=CT[b["sys"]].get(b["OL"])
    metal=b["sys"].split("-")[0]; var=b["sys"].split("-MOF")[1]
    rows.append({"System":b["sys"],"Metal":metal,"Variant":var,"Bond":f'{b["IL"]}-{b["OL"]}',
                 "d_A":b["d"],"rho_BCP":b["rho"],
                 "I_q":Idr if Idr is not None else None,
                 "O_q":Odr if Odr is not None else None,
                 "I_q_over_avg":(Idr/AVG[b["sys"]]["I"]) if (Idr and AVG[b["sys"]].get("I")) else None,
                 "O_q_over_avg":(Odr/AVG[b["sys"]]["O"]) if (Odr and AVG[b["sys"]].get("O")) else None})
nmiss=sum(1 for r in rows if r["I_q"] is None or r["O_q"] is None)
print(f"merged {len(rows)} bonds; missing-charge: {nmiss}")
# write CSV
out_csv=f"{CAMP}/IO_bond_per_atom_CT.csv"
with open(out_csv,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
# write readable txt
out_txt=f"{CAMP}/IO_bond_per_atom_CT.txt"
with open(out_txt,"w") as f:
    f.write("Per-atom charge transfer ($dr, e) on the I-O bond network\n")
    f.write("$dr = Bader-basin integral of deltarho (loaded-empty-I2only). + = atom gained e (reduced), - = lost e (oxidized).\n")
    f.write("I_q / O_q = charge on the specific I / O atom in the bond. avg given per system for reference.\n\n")
    f.write(f"{'System':14s}{'Bond':10s}{'d(Å)':>6s}{'ρ_BCP':>8s}{'I_q':>9s}{'O_q':>9s}{'I/avg':>7s}{'O/avg':>7s}  note\n")
    f.write("-"*80+"\n")
    for r in rows:
        if r["I_q"] is None or r["O_q"] is None: continue
        avgI=AVG[r["System"]]["I"]; avgO=AVG[r["System"]]["O"]
        note=""
        if abs(r["O_q"])>2*abs(avgO): note="O donates >2x avg"
        if r["I_q"]>1.5*avgI: note=(note+"; " if note else "")+"I most-charged"
        f.write(f"{r['System']:14s}{r['Bond']:10s}{r['d_A']:6.2f}{r['rho_BCP']:8.4f}{r['I_q']:+9.3f}{r['O_q']:+9.3f}{r['I_q']/avgI:7.2f}{r['O_q']/avgO:7.2f}  {note}\n")
    f.write("\nPer-system element averages (e):\n")
    for sys in SYS2DIR:
        f.write(f"  {sys:12s} avg I={AVG[sys]['I']:+.4f}  avg O={AVG[sys]['O']:+.4f}\n")
print(f"wrote {out_csv} and {out_txt}")
# ---- plot ----
fig,axes=plt.subplots(1,2,figsize=(13,5.5))
varcolor={"st":"#1f77b4","act":"#2ca02c","vc":"#d62728"}
metmark={"Y":"o","Tb":"s","Gd":"^"}
for ax,yc,ylabel,title,refkey in [
    (axes[0], lambda r:-r["O_q"], r"$-q_O$  (O donation, e)", "O donation vs I-O bond strength", "O"),
    (axes[1], lambda r:r["I_q"], r"$q_I$  (I charge, e)", "Iodine charge vs I-O bond strength", "I")]:
    for r in rows:
        if r["I_q"] is None or r["O_q"] is None: continue
        ax.scatter(r["rho_BCP"], yc(r), c=varcolor[r["Variant"]], marker=metmark[r["Metal"]],
                   s=60,alpha=0.75,edgecolor='k',linewidth=0.5,zorder=3)
    # per-variant average reference (dashed)
    for var,col in varcolor.items():
        ys=[]
        for sys in [s for s in SYS2DIR if sys.split("-MOF")[1]==var]:
            for r in rows:
                if r["System"]==sys and r["I_q"] is not None:
                    ys.append(-r["O_q"] if refkey=="O" else r["I_q"])
        if ys: ax.axhline(np.mean(ys),color=col,ls='--',lw=1,alpha=0.5,zorder=1)
    ax.set_xlabel(r"$\rho_{BCP}$ at I-O (e/bohr$^3$)",fontsize=12)
    ax.set_ylabel(ylabel,fontsize=12); ax.set_title(title,fontsize=12)
    ax.grid(alpha=0.3,zorder=0)
# legend
handles=[Line2D([0],[0],marker='o',color='w',markerfacecolor='gray',label='Y',markersize=9),
         Line2D([0],[0],marker='s',color='w',markerfacecolor='gray',label='Tb',markersize=9),
         Line2D([0],[0],marker='^',color='w',markerfacecolor='gray',label='Gd',markersize=9),
         Line2D([0],[0],marker='o',color='w',markerfacecolor=varcolor['st'],label='st',markersize=9),
         Line2D([0],[0],marker='o',color='w',markerfacecolor=varcolor['act'],label='act',markersize=9),
         Line2D([0],[0],marker='o',color='w',markerfacecolor=varcolor['vc'],label='vc',markersize=9)]
axes[1].legend(handles=handles,fontsize=9,ncol=2,loc='upper left')
fig.suptitle("Charge transfer is localized at the I-O anchoring sites (per-atom Bader $dr$)",fontsize=13)
plt.tight_layout(rect=[0,0,1,0.95])
out_png=f"{CAMP}/IO_bond_CT_vs_strength.png"
fig.savefig(out_png,dpi=150,bbox_inches='tight')
print(f"wrote {out_png}")
# ---- companion plot vs I-O distance (gradient shade + bond labels) ----
fig,axes=plt.subplots(1,2,figsize=(13,5.5))
import numpy as np
_gx=np.linspace(2.6,3.2,300)
_gv=np.where(_gx<=3.10,0.85,0.85+(_gx-3.10)/0.10*0.15)   # solid 0.85 gray to 3.1, fade to white by 3.2
_grad=_gv.reshape(1,-1)
for ax,yc,ylabel,title in [
    (axes[0], lambda r:-r["O_q"], r"$-q_O$  (O donation, e)", "O donation vs I-O distance"),
    (axes[1], lambda r:r["I_q"], r"$q_I$  (I charge, e)", "Iodine charge vs I-O distance")]:
    for r in rows:
        if r["I_q"] is None or r["O_q"] is None: continue
        ax.scatter(r["d_A"], yc(r), c=varcolor[r["Variant"]], marker=metmark[r["Metal"]],
                   s=60,alpha=0.75,edgecolor='k',linewidth=0.5,zorder=3)
    ys=[yc(r) for r in rows if r["I_q"] is not None]
    ylo,yhi=min(ys),max(ys); pad=0.08*(yhi-ylo)
    ax.imshow(_grad,aspect='auto',cmap='gray',vmin=0,vmax=1,extent=[2.6,3.2,ylo-pad,yhi+pad],
              zorder=0,origin='lower')
    ax.set_ylim(ylo-pad,yhi+pad); ax.set_xlim(2.6,4.10)
    # label each contact inside the shaded zone (d<=3.1), staggered vertically
    strong=sorted([r for r in rows if r["I_q"] and r["O_q"] and r["d_A"]<=3.10],key=lambda r:r["d_A"])
    for i,r in enumerate(strong):
        dy=[7,-10,8,-9,7,-10,8][i%7]
        if i==0:  # leftmost strong bond (I5-O12, d=2.77): label to the LEFT so it clears the next marker
            ax.annotate(r["Bond"],(r["d_A"],yc(r)),xytext=(-4,dy),textcoords='offset points',
                        fontsize=6.8,color='0.1',zorder=6,ha='right')
        else:
            ax.annotate(r["Bond"],(r["d_A"],yc(r)),xytext=(5,dy),textcoords='offset points',
                        fontsize=6.8,color='0.1',zorder=6)
    ax.set_xlabel("I-O distance (Å)",fontsize=12)
    ax.set_ylabel(ylabel,fontsize=12); ax.set_title(title,fontsize=12)
    ax.grid(alpha=0.25,zorder=1)
axes[0].legend(handles=handles,fontsize=9,ncol=2,loc='lower left',framealpha=0.9)
fig.suptitle("Per-atom charge transfer vs I-O distance (vc iodines stay highly charged even at long, weak contacts)",
             fontsize=12.5)
plt.tight_layout(rect=[0,0,1,0.95])
out_png2=f"{CAMP}/IO_bond_CT_vs_distance.png"
fig.savefig(out_png2,dpi=150,bbox_inches='tight')
print(f"wrote {out_png2}")
# ---- summary stat: correlation ----
rh=np.array([r["rho_BCP"] for r in rows if r["I_q"] and r["O_q"]])
od=np.array([-r["O_q"] for r in rows if r["I_q"] and r["O_q"]])
iq=np.array([r["I_q"] for r in rows if r["I_q"] and r["O_q"]])
print(f"\ncorr(rho_BCP, -O_q) = {np.corrcoef(rh,od)[0,1]:+.3f}")
print(f"corr(rho_BCP,  I_q) = {np.corrcoef(rh,iq)[0,1]:+.3f}")
print(f"fraction of bonds where bonded-O donates more than avg-O: {np.mean(od > np.array([abs(AVG[r['System']]['O']) for r in rows if r['I_q'] and r['O_q']])):.0%}")
