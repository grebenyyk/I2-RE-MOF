#!/usr/bin/env python3
"""Inferno-cmap variants of the publication ELF figures.
Saves files with '_inferno' suffix; does not overwrite existing magma outputs.
"""
import numpy as np, re, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import map_coordinates
plt.rcParams.update({
    "font.family":"Liberation Sans","font.size":10,
    "axes.labelsize":12,"axes.titlesize":11,
    "xtick.labelsize":10,"ytick.labelsize":10,
    "axes.linewidth":1.0,"lines.linewidth":1.8,"legend.frameon":False,
    "svg.fonttype":"none","pdf.fonttype":42,
})
BOHR=0.5291772109
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
XL=r"$e_1$ along chain, Å"; YL=r"$e_2$, Å"
CMAP='inferno'

def read_cube(fn):
    L=open(fn).readlines(); nat=int(L[2].split()[0])
    nx,ny,nz=[int(L[i].split()[0]) for i in (3,4,5)]
    ga=np.array([float(x) for x in L[3].split()[1:4]]); gb=np.array([float(x) for x in L[4].split()[1:4]]); gc=np.array([float(x) for x in L[5].split()[1:4]])
    d=np.array([float(x) for l in L[6+nat:] for x in l.split()]).reshape((nx,ny,nz))
    return d,nx,ny,nz,np.array([ga*nx,gb*ny,gc*nz]).T

def atoms(s):
    t=open(f"{BASE}/{s}/scf_loaded.in").read()
    cellA=np.array([list(map(float,re.findall(r"-?\d+\.\d+",l))) for l in t.split("CELL_PARAMETERS")[1].splitlines()[1:4]])
    apos=[]
    for l in t.split("ATOMIC_POSITIONS")[1].split("K_POINTS")[0].splitlines():
        m=re.match(r"\s*([A-Z][a-z]?)\s+(-?\d.\d+)\s+(-?\d.\d+)\s+(-?\d.\d+)",l)
        if m: apos.append((m.group(1),np.array([float(m.group(2)),float(m.group(3)),float(m.group(4))])))
    return cellA,apos

def slice_sys(s):
    elf,nx,ny,nz,cellB=read_cube(f"{BASE}/{s}/elf.cube")
    cellA,apos=atoms(s); cellT=cellA.T; invB=np.linalg.inv(cellB)
    chain=np.array([f for sp,f in apos if sp=='I' and f[2]>0.5])
    if len(chain)<3: chain=np.array([f for sp,f in apos if sp=='I'])[:4]
    cartA=(cellA@chain.T).T; cenA=cartA.mean(0)
    _,_,vh=np.linalg.svd(cartA-cenA); e1,e2,nrm=vh[0],vh[1],vh[2]; cenB=cenA/BOHR
    R=5.5; N=600; sg=np.linspace(-R,R,N)/BOHR; S,T=np.meshgrid(sg,sg)
    pts=cenB[None,None,:]+S[:,:,None]*e1[None,None,:]+T[:,:,None]*e2[None,None,:]
    fr=np.einsum('ij,kli->klj',invB,pts); fr-=np.floor(fr)
    vals=map_coordinates(elf,[fr[:,:,0].ravel()*nx,fr[:,:,1].ravel()*ny,fr[:,:,2].ravel()*nz],order=1,mode='wrap').reshape((N,N))
    atlist=[]
    for sp,f in apos:
        c=cellA@f; dd=c-cenA; atlist.append((sp,dd@e1,dd@e2,abs(dd@nrm)))
    return vals,atlist,R

def draw(ax,vals,atlist,R,vmin):
    im=ax.imshow(vals,origin='lower',extent=[-R,R,-R,R],cmap=CMAP,vmin=vmin,vmax=1.0,interpolation='bilinear')
    ax.set_aspect('equal'); ax.set_xlabel(XL); ax.set_ylabel(YL)
    for sp,ss,tt,h in atlist:
        if abs(ss)<R and abs(tt)<R and h<0.5:
            ax.scatter([ss],[tt],s=42 if sp=='I' else 14,facecolor='none',
                       edgecolor='#39FF14' if sp=='I' else('red' if sp=='O' else('cyan' if sp in('Y','Tb','Gd') else 'white')),
                       lw=1.1 if sp=='I' else 0.7,zorder=5)
            if sp=='I': ax.text(ss,tt,'I',color='white',fontsize=6.5,fontweight='bold',ha='center',va='center',zorder=6)
    return im

def plot_2panel(s):
    vals,atlist,R=slice_sys(s)
    fig,ax=plt.subplots(1,2,figsize=(7.0,3.35))
    for a,vm in zip(ax,(0.0,0.5)):
        im=draw(a,vals,atlist,R,vm)
        cb=fig.colorbar(im,ax=a,fraction=0.046,pad=0.04); cb.set_label("ELF",fontsize=11); cb.ax.tick_params(labelsize=9)
    fig.tight_layout(pad=0.4)
    fig.savefig(f"{BASE}/{s}/{s}_elf_inferno.png",dpi=600,bbox_inches='tight')
    fig.savefig(f"{BASE}/{s}/{s}_elf_inferno.pdf",bbox_inches='tight'); plt.close(fig)
    print(f"{s}: {s}_elf_inferno.png/.pdf")

def plot_str_act(metal):
    s_str=f"{metal}_str"; s_act=f"{metal}_act"
    v1,a1,R=slice_sys(s_str); v2,a2,_=slice_sys(s_act)
    fig,ax=plt.subplots(1,2,figsize=(7.0,3.35))
    for a,(v,al,lbl) in zip(ax,[(v1,a1,"Strained"),(v2,a2,"Activated")]):
        im=draw(a,v,al,R,0.5); a.set_title(lbl)
        cb=fig.colorbar(im,ax=a,fraction=0.046,pad=0.04); cb.set_label("ELF",fontsize=11); cb.ax.tick_params(labelsize=9)
    fig.tight_layout(pad=0.4)
    fig.savefig(f"{BASE}/{metal}_elf_str_act_inferno.png",dpi=600,bbox_inches='tight')
    fig.savefig(f"{BASE}/{metal}_elf_str_act_inferno.pdf",bbox_inches='tight'); plt.close(fig)
    print(f"{metal}: {metal}_elf_str_act_inferno.png/.pdf (str|act, ELF 0.5-1)")

for s in ["Y_str","Tb_str","Gd_str","Y_act","Tb_act","Gd_act"]:
    plot_2panel(s)
for m in ["Y","Tb","Gd"]:
    plot_str_act(m)
