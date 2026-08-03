#!/usr/bin/env python3
"""4x4x2 redo for the binding tables. Fully self-consistent at 50/400 PBE+D3 mv0.01:
  loaded (F+I) and isubl (I-sublattice) @ 4x4x2  [metallic, k-redo]
  empty framework @ Gamma                         [semiconducting, k-exact]
  gasI2 @ Gamma (have), vc-relaxed empty @ Gamma (have, from vc_empty run).
Reads geometries from qtaim/campaign/<variant>/scf_{loaded,i2only,empty}.in."""
import os,re,numpy as np
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
OUT=f"{BASE}/convergence_test_Yvc"
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
     "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF",
     "Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF","Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF","N":"N.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"I":126.904,"Gd":157.25,"Tb":158.925,"N":14.007}
LOADED=["Y_str","Y_actc","Y_deep","Gd_str","Gd_act","Gd_deep","Tb_str","Tb_act","Tb_deep"]
ISUBL=["Y_str","Y_actc","Y_deep","Gd_deep","Tb_deep"]
def geom(v,kind):
    f=f"{BASE}/{v}/scf_{kind}.in"
    L=open(f).read().splitlines()
    ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
    ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][0]
    cell=np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)])
    pos=[];k=ai+1
    while k<len(L):
        t=L[k].split()
        if not t or t[0] in ("K_POINTS","ATOMIC_SPECIES"): break
        if t[0] in PSE and len(t)>=4:
            try: pos.append((t[0],float(t[1]),float(t[2]),float(t[3])))
            except: pass
        k+=1
    return cell,pos
def emit(tag,cell,pos,kmesh):
    uniq=[]
    for e in [p[0] for p in pos]:
        if e not in uniq: uniq.append(e)
    prefix=f"redo442_{tag}"
    L=[f"&CONTROL\n  calculation = 'scf'\n  prefix = '{prefix}'\n  outdir = './tmp_{prefix}'",
       "  pseudo_dir = '/workspace/qe-7.5/pseudo'\n  tprnfor = .true.\n/",
       f"&SYSTEM\n  ibrav = 0\n  nat = {len(pos)}\n  ntyp = {len(uniq)}\n  ecutwfc = 50\n  ecutrho = 400",
       "  occupations = 'smearing'\n  smearing = 'mv'\n  degauss = 0.01\n  nspin = 1\n  vdw_corr = 'grimme-d3'\n  dftd3_version = 4\n/",
       "&ELECTRONS\n  conv_thr = 1.0d-8\n  diagonalization = 'david'\n  diago_david_ndim = 2",
       "  mixing_mode = 'local-TF'\n  mixing_beta = 0.2\n  mixing_ndim = 16\n  electron_maxstep = 300\n/",
       "ATOMIC_SPECIES"]
    for e in uniq: L.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell: L.append(f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for e,x,y,z in pos: L.append(f"  {e:<2s} {x:.10f} {y:.10f} {z:.10f}")
    L.append(f"K_POINTS automatic\n  {kmesh}\n")
    open(f"{OUT}/{prefix}.in","w").write("\n".join(L)); return prefix
man=[]
for v in LOADED:
    c,p=geom(v,"loaded"); man.append(emit(f"loaded_{v}",c,p,"4 4 2 0 0 0"))
for v in ISUBL:
    c,p=geom(v,"i2only"); man.append(emit(f"isubl_{v}",c,p,"4 4 2 0 0 0"))
for v in LOADED:
    c,p=geom(v,"empty"); man.append(emit(f"empty_{v}",c,p,"1 1 1 0 0 0"))
open(f"{OUT}/redo442_manifest.txt","w").write("\n".join(man)+"\n")
print(f"wrote {len(man)} inputs:")
for t in man: print("  ",t)
