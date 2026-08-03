import numpy as np, re
from collections import defaultdict
CUBE="13Y_diffdens.cube"; SCF="13Y_opt/scf.in"; BOHR=0.529177
RAD={"H":1.20,"C":1.70,"O":1.52,"Y":2.00,"I":1.98}  # vdW Å
lines=open(CUBE).readlines(); nat=int(lines[2].split()[0])
nx,ny,nz=[int(lines[i].split()[0]) for i in (3,4,5)]
ga=np.array([float(x) for x in lines[3].split()[1:4]]); gb=np.array([float(x) for x in lines[4].split()[1:4]]); gc=np.array([float(x) for x in lines[5].split()[1:4]])
a=ga*nx; b=gb*ny; c=gc*nz; cell=np.array([a,b,c])
data=np.array([float(x) for l in lines[6+nat:] for x in l.split()]).reshape((nx,ny,nz))
vox=abs(np.dot(ga,np.cross(gb,gc)))
s=open(SCF).read()
apos=[]
for l in s.split("ATOMIC_POSITIONS")[1].split("K_POINTS")[0].splitlines():
    m=re.match(r"\s*([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",l)
    if m: apos.append((m.group(1),np.array([float(m.group(2)),float(m.group(3)),float(m.group(4))])))
FX,FY,FZ=np.meshgrid(np.arange(nx),np.arange(ny),np.arange(nz),indexing='ij')
q=defaultdict(float); qcount=defaultdict(int)
for sym,fr in apos:
    dfx=FX/nx-fr[0]; dfy=FY/ny-fr[1]; dfz=FZ/nz-fr[2]
    dfx-=np.round(dfx); dfy-=np.round(dfy); dfz-=np.round(dfz)
    dcart=dfx[...,None]*a+dfy[...,None]*b+dfz[...,None]*c
    r=np.linalg.norm(dcart,axis=-1); m=r<RAD[sym]/BOHR
    q[sym]+=data[m].sum()*vox; qcount[sym]+=1
print(f"vdW-sphere integrated δρ (per species, summed):")
print(f"{'spec':>5}{'N':>4}{'Σδρ(e)':>10}{'per-atom':>10}")
tot=0
for s_ in ['C','H','O','Y','I']:
    n=qcount[s_]; per=q[s_]/n if n else 0; print(f"{s_:>5}{n:>4}{q[s_]:>+10.3f}{per:>+10.3f}"); tot+=q[s_]
print(f"{'Σ spheres':>9}{tot:>+10.3f}  (residual = CT outside vdW spheres: polariz./bond density)")
print(f"\nI gained {q['I']:+.2f} e ; donors: O {q['O']:+.2f}, C {q['C']:+.2f}, Y {q['Y']:+.2f}, H {q['H']:+.2f}")
print(f"donor ranking: {' > '.join(f'{s}({q[s]:+.2f})' for s in sorted(['C','H','O','Y'],key=lambda x:q[x]))}")
