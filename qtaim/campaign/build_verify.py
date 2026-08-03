#!/usr/bin/env python3
"""Verification suite (Aspect 1 vc-clean + Aspect 2 MP-smearing + Aspect 3 molecular-vs-polyiodide).
All on the EXACT production geometries (scf_loaded.in / scf_i2only.in / activated_i2_mol relax-final).
  vc loaded+isubl @ 60/800 4x4x2  (clean vc k-test vs the 60/800 Gamma scf_loaded/i2only)
  MP-smearing @ 4x4x2 (Y-act 50/400, Y-vc 60/800)  (Aspect 2: rule out mv artifact)
  molecular-I2 @ 50/400 4x4x2 (Y-act)  (Aspect 3: polyiodide vs molecular at proper k)"""
import os,re,numpy as np
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"; ROOT="/workspace/dimitry/calc/qe/i2_mof"
OUT=f"{BASE}/convergence_test_Yvc"
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF","Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF","Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF","Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF","N":"N.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"I":126.904,"Gd":157.25,"Tb":158.925,"N":14.007}
def geom_inp(f):
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
def geom_out_final(f):
    t=open(f).read()
    m=re.search(r"CELL_PARAMETERS \(angstrom\)\s*\n((?:\s*-?\d[\d.\s-]*\n){3})",t)
    cell=np.array([[float(x) for x in l.split()[:3]] for l in m.group(1).splitlines()]) if m else None
    m2=re.search(r"Begin final coordinates.*?ATOMIC_POSITIONS\s*\(crystal\)\s*\n(.*?)End final coordinates",t,re.DOTALL)
    pos=[]
    for l in m2.group(1).splitlines():
        s=l.split()
        if len(s)>=4 and s[0] in PSE:
            try: pos.append((s[0],float(s[1]),float(s[2]),float(s[3])))
            except: pass
    return cell,pos
def emit(tag,cell,pos,ew,er,kmesh,smear):
    uniq=[]
    for e in [p[0] for p in pos]:
        if e not in uniq: uniq.append(e)
    mb=0.1 if smear=='mp' else 0.2
    L=[f"&CONTROL\n  calculation = 'scf'\n  prefix = '{tag}'\n  outdir = './tmp_{tag}'",
       "  pseudo_dir = '/workspace/qe-7.5/pseudo'\n  tprnfor = .true.\n/",
       f"&SYSTEM\n  ibrav = 0\n  nat = {len(pos)}\n  ntyp = {len(uniq)}\n  ecutwfc = {ew}\n  ecutrho = {er}",
       f"  occupations = 'smearing'\n  smearing = '{smear}'\n  degauss = 0.01\n  nspin = 1\n  vdw_corr = 'grimme-d3'\n  dftd3_version = 4\n/",
       "&ELECTRONS\n  conv_thr = 1.0d-8\n  diagonalization = 'david'\n  diago_david_ndim = 2",
       f"  mixing_mode = 'local-TF'\n  mixing_beta = {mb}\n  mixing_ndim = 24\n  electron_maxstep = 400\n/",
       "ATOMIC_SPECIES"]
    for e in uniq: L.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell: L.append(f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for e,x,y,z in pos: L.append(f"  {e:<2s} {x:.10f} {y:.10f} {z:.10f}")
    L.append(f"K_POINTS automatic\n  {kmesh}\n")
    open(f"{OUT}/{tag}.in","w").write("\n".join(L)); return tag
man=[]
# 1) vc loaded + isubl @ 60/800 4x4x2 (clean vc k-test)
for v in ["Y_deep","Gd_deep","Tb_deep"]:
    c,p=geom_inp(f"{BASE}/{v}/scf_loaded.in"); man.append(emit(f"verify_vcL_{v}",c,p,60,800,"4 4 2 0 0 0","mv"))
    c,p=geom_inp(f"{BASE}/{v}/scf_i2only.in"); man.append(emit(f"verify_vcI_{v}",c,p,60,800,"4 4 2 0 0 0","mv"))
# 2) MP-smearing spot-checks
c,p=geom_inp(f"{BASE}/Y_actc/scf_loaded.in"); man.append(emit("verify_mp_yact",c,p,50,400,"4 4 2 0 0 0","mp"))
c,p=geom_inp(f"{BASE}/Y_deep/scf_loaded.in"); man.append(emit("verify_mp_yvc",c,p,60,800,"4 4 2 0 0 0","mp"))
# 3) molecular-I2 (Aspect 3): Y-act activated_i2_mol relax-final @ 50/400 4x4x2
molf=f"{ROOT}/13_Y/activated_i2_mol/loaded.out"
if os.path.exists(molf):
    c,p=geom_out_final(molf); man.append(emit("verify_mol_yact",c,p,50,400,"4 4 2 0 0 0","mv"))
else: print(f"[warn] molecular geom not found: {molf}")
open(f"{OUT}/verify_manifest.txt","w").write("\n".join(man)+"\n")
print(f"wrote {len(man)} verify inputs:")
for t in man: print("  ",t)
