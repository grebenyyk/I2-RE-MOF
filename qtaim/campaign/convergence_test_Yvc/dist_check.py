#!/usr/bin/env python3
"""I-I min-image contact check for the HSE poly/mol geometries (act, 3 metals).
Molecular reference should show 4 isolated ~2.7A pairs (8 I = 4 I2), nothing else <3.5A.
Polyiodide should show a chain (~7 contacts ~2.7-3.0A). Flags a bad molecular reference."""
import sys,numpy as np
def read_geom(fn):
    L=open(fn).read().splitlines()
    ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
    cell=np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)])
    ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][0]
    fr=[];k=ai+1
    while k<len(L):
        t=L[k].split()
        if not t or t[0] in ("K_POINTS","ATOMIC_SPECIES"): break
        if t[0]=="I": fr.append([float(t[1]),float(t[2]),float(t[3])])
        k+=1
    return cell,np.array(fr)
def contacts(cell,fr,cut=4.0):
    out=[]
    for i in range(len(fr)):
        for j in range(i+1,len(fr)):
            df=fr[i]-fr[j]; df-=np.round(df)
            d=np.linalg.norm(df@cell)
            if d<cut: out.append((i,j,d))
    return sorted(out,key=lambda x:x[2])
for fn in sys.argv[1:]:
    cell,fr=read_geom(fn)
    c=contacts(cell,fr)
    short=[d for *_,d in c if d<3.25]
    print(f"\n{fn}: {len(fr)} I, {len(c)} contacts<4.0A, {len(short)} short(<3.25):")
    for i,j,d in c: print(f"  I{i}-I{j}: {d:.3f}")
