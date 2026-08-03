#!/usr/bin/env python3
"""13Y ELF — two-panel slice through the I4 polyiodide plane.
Left: ELF 0->1. Right: ELF>=0.5 (localized: bonds/lone pairs).
Atoms overlaid: lime=I (labeled), red=O, white=C/H, cyan=Y."""
import numpy as np, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import map_coordinates
CUBE="elf.cube"; SCF="scf.in"; BOHR=0.529177
L=open(CUBE).readlines(); nat=int(L[2].split()[0])
nx,ny,nz=[int(L[i].split()[0]) for i in (3,4,5)]
ga=np.array([float(x) for x in L[3].split()[1:4]]); gb=np.array([float(x) for x in L[4].split()[1:4]]); gc=np.array([float(x) for x in L[5].split()[1:4]])
cellB=np.array([ga*nx,gb*ny,gc*nz]).T; invB=np.linalg.inv(cellB)
data=np.array([float(x) for l in L[6+nat:] for x in l.split()]).reshape((nx,ny,nz))
apos=[]
for l in open(SCF).read().split("ATOMIC_POSITIONS")[1].split("K_POINTS")[0].splitlines():
    m=re.match(r"\s*([A-Z][a-z]?)\s+(-?\d.\d+)\s+(-?\d.\d+)\s+(-?\d.\d+)",l)
    if m: apos.append((m.group(1),np.array([float(m.group(2)),float(m.group(3)),float(m.group(4))])))
cellA=np.array([[8.9348,0,0],[0,9.6608,0],[-3.76256,0,16.59583]])
chain=np.array([f for s,f in apos if s=='I' and f[2]>0.5]); cartA=(cellA@chain.T).T; cenA=cartA.mean(0)
_,_,vh=np.linalg.svd(cartA-cenA); e1=vh[0]; e2=vh[1]; nrm=vh[2]; cenB=cenA/BOHR
R=5.5; N=700; sg=np.linspace(-R,R,N); sB=sg/BOHR; S,T=np.meshgrid(sB,sB)
pts=cenB[None,None,:]+S[:,:,None]*e1[None,None,:]+T[:,:,None]*e2[None,None,:]
fr=np.einsum('ij,kli->klj',invB,pts); fr-=np.floor(fr)
vals=map_coordinates(data,[fr[:,:,0].ravel()*nx,fr[:,:,1].ravel()*ny,fr[:,:,2].ravel()*nz],order=1,mode='wrap').reshape((N,N))
atlist=[]
for sym,f in apos:
    c=cellA@f; d=c-cenA; atlist.append((sym,d@e1,d@e2,abs(d@nrm)))
fig,ax=plt.subplots(1,2,figsize=(12,5.8))
for a,(vm) in zip(ax,(0.0,0.5)):
    im=a.imshow(vals,origin='lower',extent=[-R,R,-R,R],cmap='magma',vmin=vm,vmax=1.0,interpolation='bilinear')
    a.set_aspect('equal'); a.set_xlabel('e1 along chain (Å)'); a.set_ylabel('e2 (Å)')
    for sym,s,t,h in atlist:
        if abs(s)<R and abs(t)<R and h<0.5:
            a.scatter([s],[t],s=70 if sym=='I' else 25,facecolor='none',
                      edgecolor='lime' if sym=='I' else('red' if sym=='O' else 'white'),
                      lw=1.6 if sym=='I' else 1.0,zorder=5)
            if sym=='I': a.text(s,t,'  I',color='white',fontsize=8,fontweight='bold',va='center')
    plt.colorbar(im,ax=a,fraction=0.046,label='ELF')
plt.tight_layout()
plt.savefig("13Y_elf_chain_plane.png",dpi=160,bbox_inches='tight')
print("wrote 13Y_elf_chain_plane.png (two panels, no title/profile)")
