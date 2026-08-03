#!/usr/bin/env python3
"""Per-atom CT on ALL geometric I-O contacts (not just BCP-having bonds).
Self-consistent: distances AND charges both measured on the QTIM cube (rho17.cube),
which is the structure the Bader-basin deltarho charges live on."""
import re,os,numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
CAMP="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"; BOHR=0.529177
SYS2DIR={"Y-MOFst":"Y_str","Tb-MOFst":"Tb_str","Gd-MOFst":"Gd_str",
         "Y-MOFact":"Y_actc","Tb-MOFact":"Tb_act","Gd-MOFact":"Gd_act",
         "Y-MOFvc":"Y_deep","Tb-MOFvc":"Tb_deep","Gd-MOFvc":"Gd_deep"}
def parse_cube(fn):
    L=open(fn,errors='ignore').read().splitlines()
    natom=int(L[2].split()[0]); bohr=natom>0; natom=abs(natom)
    cell=np.zeros((3,3))
    for k in range(3):
        q=L[3+k].split(); n=int(q[0]); cell[k]=np.array([float(q[1]),float(q[2]),float(q[3])])*n
    if bohr: cell=cell*BOHR
    atoms=[]
    for k in range(natom):
        q=L[6+k].split(); Z=int(float(q[0])); xyz=np.array([float(q[2]),float(q[3]),float(q[4])])
        if bohr: xyz=xyz*BOHR
        atoms.append((Z,xyz))
    return cell,atoms
def parse_ct(d):
    rows=[]
    pat=re.compile(r'^\s*(\d+)\s+\d+\s+\d+\s+([A-Z][a-z]?)\s+\d+\s+--\s+([0-9.+\-E]+)\s+([0-9.+\-E]+)\s+([0-9.+\-E]+)\s+([0-9.+\-E]+)')
    for l in open(f"{CAMP}/{d}/ct.out",errors="ignore"):
        m=pat.match(l)
        if m: rows.append((int(m.group(1)),m.group(2),float(m.group(6))))
    rows.sort(key=lambda r:r[0]); lbl={}; cnt={}
    for iat,el,dr in rows:
        cnt[el]=cnt.get(el,0)+1; lbl[f"{el}{cnt[el]}"]=dr
    return lbl
CUTOFF=4.0
recs=[]
print(f"{'system':10s}{'shortest':>9s}{'#<3.1':>7s}{'#<4.0':>7s}{'contacts':>10s}")
for sys,d in SYS2DIR.items():
    cell,atoms=parse_cube(f"{CAMP}/{d}/rho17.cube")
    inv=np.linalg.inv(cell)
    # within-element labels by cube atom order (matches ct.out basin order)
    cnt={}; labeled=[]
    for Z,xyz in atoms:
        el={6:'C',1:'H',8:'O',39:'Y',65:'Tb',64:'Gd',53:'I'}.get(Z,'X')
        cnt[el]=cnt.get(el,0)+1; labeled.append((f"{el}{cnt[el]}",xyz))
    ct=parse_ct(d)
    I=[(lab,xyz) for lab,xyz in labeled if lab.startswith('I')]
    O=[(lab,xyz) for lab,xyz in labeled if lab.startswith('O')]
    metal=sys.split('-')[0]; var=sys.split('-MOF')[1]; ds_all=[]; n=0
    for ilab,ip in I:
        for olab,op in O:
            delta=op-ip; frac=delta@inv; frac=frac-np.round(frac); dist=np.linalg.norm(frac@cell)
            ds_all.append(dist)
            if dist<CUTOFF and ilab in ct and olab in ct:
                recs.append({"Metal":metal,"Variant":var,"d":dist,"I_q":ct[ilab],"O_q":ct[olab],"Bond":f"{ilab}-{olab}"})
                n+=1
    ds_all=np.array(sorted(ds_all))
    print(f"{sys:10s}{ds_all[0]:7.2f} Å {np.sum(ds_all<3.1):7d}{np.sum(ds_all<4.0):7d}{n:10d}")
print(f"\ntotal geometric I-O contacts within {CUTOFF} Å: {len(recs)}")
# ---- plot ----
varcolor={"st":"#1f77b4","act":"#2ca02c","vc":"#d62728"}; metmark={"Y":"o","Tb":"s","Gd":"^"}
fig,axes=plt.subplots(1,2,figsize=(13,5.5))
_gx=np.linspace(2.6,3.2,300)
_gv=np.where(_gx<=3.10,0.85,0.85+(_gx-3.10)/0.10*0.15)
_grad=_gv.reshape(1,-1)
for ax,yc,ylabel,title in [
    (axes[0], lambda r:-r["O_q"], r"$-q_O$  (O donation, e)", "O donation vs I-O distance"),
    (axes[1], lambda r:r["I_q"], r"$q_I$  (I charge, e)", "Iodine charge vs I-O distance")]:
    for r in recs:
        ax.scatter(r["d"], yc(r), c=varcolor[r["Variant"]], marker=metmark[r["Metal"]],
                   s=35,alpha=0.6,edgecolor='k',linewidth=0.3,zorder=3)
    ys=[yc(r) for r in recs]; ylo,yhi=min(ys),max(ys); pad=0.08*(yhi-ylo)
    ax.imshow(_grad,aspect='auto',cmap='gray',vmin=0,vmax=1,extent=[2.6,3.2,ylo-pad,yhi+pad],zorder=0,origin='lower')
    ax.set_ylim(ylo-pad,yhi+pad); ax.set_xlim(2.6,4.05)
    ax.set_xlabel("I-O distance (Å)",fontsize=12); ax.set_ylabel(ylabel,fontsize=12); ax.set_title(title,fontsize=12)
    ax.grid(alpha=0.25,zorder=1)
handles=[Line2D([0],[0],marker='o',color='w',markerfacecolor='gray',label='Y',markersize=9),
         Line2D([0],[0],marker='s',color='w',markerfacecolor='gray',label='Tb',markersize=9),
         Line2D([0],[0],marker='^',color='w',markerfacecolor='gray',label='Gd',markersize=9),
         Line2D([0],[0],marker='o',color='w',markerfacecolor=varcolor['st'],label='st',markersize=9),
         Line2D([0],[0],marker='o',color='w',markerfacecolor=varcolor['act'],label='act',markersize=9),
         Line2D([0],[0],marker='o',color='w',markerfacecolor=varcolor['vc'],label='vc',markersize=9)]
axes[0].legend(handles=handles,fontsize=9,ncol=2,loc='lower left',framealpha=0.9)
fig.suptitle("Per-atom charge transfer vs I-O distance (all geometric contacts; cube=self-consistent)",fontsize=12.5)
plt.tight_layout(rect=[0,0,1,0.95])
out=f"{CAMP}/IO_bond_CT_vs_distance_geometric.png"
fig.savefig(out,dpi=150,bbox_inches='tight'); print(f"wrote {out}")
rh=np.array([r["d"] for r in recs]); od=np.array([-r["O_q"] for r in recs]); iq=np.array([r["I_q"] for r in recs])
print(f"corr(d, O donation)={np.corrcoef(rh,od)[0,1]:+.3f}  corr(d, I charge)={np.corrcoef(rh,iq)[0,1]:+.3f}")
