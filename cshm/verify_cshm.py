#!/usr/bin/env python3
"""Verify the CShM table: recompute CN=6 CShM for vc (unverified) + activated (spot-check),
and CN=8 for as-synthesized, using the same protocol as cshm_i2mof.py (CUT_O=2.75, cosymlib)."""
import warnings, numpy as np; warnings.filterwarnings("ignore")
from ase.io import read
from ase.neighborlist import neighbor_list
import spglib
from cosymlib import Geometry
from cosymlib.shape import Shape
from cosymlib.shape.tools import get_structure_references
CUT=2.75; DONORS={"O","N","S","P","F","Cl","Br"}; METALS=set("Y Gd Tb La Ce Pr Nd Sm Eu Dy Ho Er Tm Yb Lu Sc Zr Hf Ti Zn".split())
FR="/workspace/dimitry/calc/qe/i2_mof/final_activated_results"
CSD="/workspace/dimitry/calc/qe/i2_mof/cshm"
def donors(at,idx,cut):
    iA,jA,dA,vA=neighbor_list("ijdD",at,cut+0.05); out=[]
    for d,j,v in zip(dA[iA==idx],jA[iA==idx],vA[iA==idx]):
        if j!=idx and at[j].symbol in DONORS and d<=cut: out.append((float(d),int(j),np.asarray(v)))
    out.sort(key=lambda x:x[0]); return out
def uniq(at,midx):
    ds=spglib.get_symmetry_dataset((np.asarray(at.get_cell()),at.get_scaled_positions(),np.asarray(at.numbers)),symprec=1e-2)
    eq=ds.equivalent_atoms; seen={}; reps=[]
    for mi in midx:
        k=int(eq[mi])
        if k not in seen: seen[k]=mi; reps.append(mi)
    return reps
def topn(coords,cn,n=5):
    refs=[str(r) for r in get_structure_references(cn)]
    c=coords-coords.mean(axis=0); shp=Shape(Geometry(c,["X"]*cn)); rows=[]
    for lab in refs:
        try: rows.append((float(shp.measure(lab,central_atom=0)),lab))
        except: pass
    rows.sort(); return rows[:n]
def run(name,at,cn):
    at.pbc=True; midx=[i for i,s in enumerate(at.get_chemical_symbols()) if s in METALS]
    reps=uniq(at,midx); pos=at.get_positions()
    for mi in reps:
        dn=donors(at,mi,CUT); sym=at[mi].symbol
        if len(dn)!=cn: print(f"  {name} {sym}: CN={len(dn)} (expected {cn}) [SKIP]"); continue
        coords=np.array([pos[mi]+dn[k][2] for k in range(len(dn))])
        rows=topn(coords,cn)
        print(f"  {name} {sym} CN={cn}: "+ "  ".join(f"{lab}({s:.4f})" for s,lab in rows))

print("=== ACTIVATED (re-confirm) ===")
for el in ["Gd","Tb","Y"]:
    try: run(f"act-{el}", read(f"{CSD}/{el}_POSCAR_activated_empty"), 6)
    except Exception as e: print(f"  act-{el}: {e}")
print("\n=== VARIABLE-CELL (verify) ===")
for el in ["Gd","Tb","Y"]:
    try: run(f"vc-{el}", read(f"{FR}/vc-relaxed/{el}_vc_empty.cif"), 6)
    except Exception as e: print(f"  vc-{el}: {e}")
print("\n=== AS-SYNTHESIZED (re-confirm, CN=8) ===")
for at in read(f"{CSD}/1907532-1999761.cif", index=":"):
    el=[s for s in at.get_chemical_symbols() if s in METALS]
    if el: run(f"pristine-{el[0]}", at, 8)
