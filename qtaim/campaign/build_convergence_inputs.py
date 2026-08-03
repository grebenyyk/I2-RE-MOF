#!/usr/bin/env python3
"""Convergence-test inputs for I2@Y-MOFvc (v2: wf_collect moved to &CONTROL;
gasI2 box 20->15 A to avoid GPU OOM at ecutrho>=800). See v1 header for full doc."""
import os, re, numpy as np

BASE = "/workspace/dimitry/calc/qe/i2_mof"
OUT  = f"{BASE}/qtaim/campaign/convergence_test_Yvc"
os.makedirs(OUT, exist_ok=True)

PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
     "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF",
     "Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF","Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF",
     "N":"N.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"I":126.904,"Gd":157.25,"Tb":158.925,"N":14.007}

COMP={
 "loaded":   (f"{BASE}/qtaim/campaign/Y_deep/scf_loaded.in", False),
 "emptyFr":  (f"{BASE}/qtaim/campaign/Y_deep/scf_empty.in",  False),
 "isubl":    (f"{BASE}/qtaim/campaign/Y_deep/scf_i2only.in", False),
 "gasI2":    (f"{BASE}/i2_vac_ref/scf.in", False),
 "emptyRel": (f"{BASE}/qtaim/campaign/Y_deep/vc_empty.out",  True),
}

def parse_geom(path, is_output):
    L=open(path).read().splitlines()
    if is_output:
        ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][-1]
        ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][-1]
    else:
        ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
        ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][0]
    cell=np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)])
    pos=[]; k=ai+1
    while k<len(L):
        t=L[k].split()
        if not t or t[0] in ("K_POINTS","ATOMIC_SPECIES"): break
        if t[0] in PSE and len(t)>=4:
            try: pos.append((t[0],float(t[1]),float(t[2]),float(t[3])))
            except ValueError: pass
        k+=1
        if len(pos)>=400: break
    return cell, pos

def emit(comp, cell, pos, ew, er, kmesh, dg, tag, pdos=False):
    elems=[p[0] for p in pos]; nat=len(pos)
    uniq=[]; [uniq.append(e) for e in elems if e not in uniq]
    ntyp=len(uniq); prefix=f"conv_{tag}"; outdir=f"./tmp_{tag}"
    mb=0.1 if dg<=0.005 else 0.2; emstep=500 if dg<=0.005 else 300
    wfc="  wf_collect = .true." if pdos else ""    # CORRECT: wf_collect is a &CONTROL var in QE7.5
    L=[f"&CONTROL\n  calculation = 'scf'\n  prefix = '{prefix}'\n  outdir = '{outdir}'",
       f"  pseudo_dir = '/workspace/qe-7.5/pseudo'\n  tprnfor = .true.{wfc}\n/",
       f"&SYSTEM\n  ibrav = 0\n  nat = {nat}\n  ntyp = {ntyp}",
       f"  ecutwfc = {ew}\n  ecutrho = {er}\n  occupations = 'smearing'",
       f"  smearing = 'mv'\n  degauss = {dg}\n  nspin = 1\n  vdw_corr = 'grimme-d3'",
       "  dftd3_version = 4\n/",
       "&ELECTRONS\n  conv_thr = 1.0d-8\n  diagonalization = 'david'",
       "  diago_david_ndim = 2\n  mixing_mode = 'local-TF'",
       f"  mixing_beta = {mb}\n  mixing_ndim = 16\n  electron_maxstep = {emstep}\n/",
       "ATOMIC_SPECIES"]
    for e in uniq: L.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell: L.append(f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for e,x,y,z in pos: L.append(f"  {e:<2s} {x:.10f} {y:.10f} {z:.10f}")
    L.append(f"K_POINTS automatic\n  {kmesh}\n")
    f=f"{OUT}/{tag}.in"; open(f,"w").write("\n".join(L)); return f

GEO={c:parse_geom(*v) for c,v in COMP.items()}
# gasI2: 20 A box OOMs the GPU at ecutrho>=800 (56 GB density). 15 A isolates I2 fine
# (nonpolar; image ~12 A away) and fits; fractional positions are box-independent for cubic.
_c,_p=GEO["gasI2"]; s=20.0/15.0   # keep Cartesian bond: frac coords scale with box (20->15 A)
GEO["gasI2"]=(15.0*np.eye(3), [(e,x*s,y*s,z*s) for (e,x,y,z) in _p])

manifest=[]
def add(comp,ew,er,kmesh,dg,tag,pdos=False):
    cell,pos=GEO[comp]; f=emit(comp,cell,pos,ew,er,kmesh,dg,tag,pdos); manifest.append((tag,f,pdos))

for c in ["loaded","emptyFr","isubl","gasI2","emptyRel"]:
    for ew,er,t in [(50,400,"5040"),(60,800,"6080"),(70,1000,"70100")]:
        add(c,ew,er,"1 1 1 0 0 0",0.01,f"{c}_{t}_Gamma")
for c in ["loaded","isubl"]:
    for kn,km in [("k211","2 1 1 0 0 0"),("k221","2 2 1 0 0 0")]:
        add(c,50,400,km,0.01,f"{c}_5040_{kn}")
for kn,km in [("k211","2 1 1 0 0 0"),("k221","2 2 1 0 0 0")]:
    add("loaded",60,800,km,0.01,f"loaded_6080_{kn}",pdos=True)
for dg in [0.005,0.002,0.001]:
    add("loaded",60,800,"2 2 1 0 0 0",dg,f"loaded_6080_k221_dg{str(dg).replace('.','p')}")

with open(f"{OUT}/manifest.txt","w") as m:
    m.write("# tag\tinput\tpdos_after?\n")
    for tag,f,pd in manifest: m.write(f"{tag}\t{f}\t{'PDOS' if pd else ''}\n")
print(f"v2: wrote {len(manifest)} inputs (wf_collect in &CONTROL; gasI2 15 A box)")
