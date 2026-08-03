#!/usr/bin/env python3
"""vc column as 4-site aggregate (P1 -> 4 inequivalent RE). For each CN6 reference
shape, compute mean and [min,max] of the CShM across the 4 RE sites; rank by mean."""
import warnings, numpy as np; warnings.filterwarnings("ignore")
from ase.io import read; from ase.neighborlist import neighbor_list
from cosymlib import Geometry; from cosymlib.shape import Shape
from cosymlib.shape.tools import get_structure_references
FR="/workspace/dimitry/calc/qe/i2_mof/final_activated_results"; CUT=2.75
def donors(at,idx):
    iA,jA,dA,vA=neighbor_list("ijdD",at,CUT+0.05); o=[]
    for d,j,v in zip(dA[iA==idx],jA[iA==idx],vA[iA==idx]):
        if j!=idx and at[j].symbol in("O","N") and d<=CUT: o.append(np.asarray(v))
    return o
print("VARIABLE-CELL (vc_empty, CN=6, 4 inequivalent RE sites — P1)\n")
print("For each metal: top-5 shapes ranked by 4-site MEAN CShM, shown as mean [min-max]\n")
for el in ["Gd","Tb","Y"]:
    at=read(f"{FR}/vc-relaxed/{el}_vc_empty.cif"); at.pbc=True; pos=at.get_positions()
    midx=[i for i,a in enumerate(at) if a.symbol in("Y","Gd","Tb")]
    refs=[str(r) for r in get_structure_references(6)]
    per_site=[]   # list of {label: S} per RE site
    for mi in midx:
        dn=donors(at,mi)
        if len(dn)!=6: print(f"  {el} atom{mi}: CN={len(dn)} (skip)"); continue
        co=np.array([pos[mi]+dn[k] for k in range(6)])-pos[mi]  # metal-centered
        co=co-co.mean(0); shp=Shape(Geometry(co,["X"]*6)); m={}
        for lab in refs:
            try: m[lab]=float(shp.measure(lab,central_atom=0))
            except: pass
        per_site.append(m)
    # aggregate per shape
    agg=[]
    for lab in refs:
        vals=[s[lab] for s in per_site if lab in s]
        if len(vals)==len(per_site): agg.append((np.mean(vals),min(vals),max(vals),lab))
    agg.sort()
    print(f"=== {el}-MOFvc ({len(per_site)} sites) ===")
    for mean,mn,mx,lab in agg[:5]:
        print(f"  {lab:8s} mean={mean:6.2f}  [{mn:.2f}-{mx:.2f}]")
    print()
