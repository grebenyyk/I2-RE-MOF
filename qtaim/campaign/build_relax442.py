#!/usr/bin/env python3
"""Full 4x4x2 relaxation of iodine only (framework frozen) from TWO starting points:
  molecular I2 (RASPA-seeded, activated_i2_mol relax-final)
  polyiodide   (activated_i2 relax-final = scf_loaded.in Y_actc)
Same activated cell, same setup (50/400, PBE+D3, mv 0.01, davidson, 4x4x2).
Framework atoms (1-72) frozen via if_pos 0 0 0; iodine (73-80) free (1 1 1).
Answers: does polyiodide form at proper k, or does it dissociate to molecular?"""
import re,os,numpy as np
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"; ROOT="/workspace/dimitry/calc/qe/i2_mof"
OUT=f"{BASE}/convergence_test_Yvc"
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF","Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"I":126.904}
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
    """Cell from the .in, positions from .out 'Begin final coordinates'."""
    inL=open(f.replace(".out",".in")).read().splitlines()
    ci=[i for i,l in enumerate(inL) if l.startswith("CELL_PARAMETERS")][0]
    cell=np.array([[float(x) for x in inL[ci+1+k].split()[:3]] for k in range(3)])
    t=open(f).read()
    m=re.search(r"Begin final coordinates.*?ATOMIC_POSITIONS\s*\(crystal\)\s*\n(.*?)End final coordinates",t,re.DOTALL)
    pos=[]
    for l in m.group(1).splitlines():
        s=l.split()
        if len(s)>=4 and s[0] in PSE:
            try: pos.append((s[0],float(s[1]),float(s[2]),float(s[3])))
            except: pass
    return cell,pos
def emit(tag,cell,pos):
    uniq=[]
    for e in [p[0] for p in pos]:
        if e not in uniq: uniq.append(e)
    nat=len(pos); nfw=nat-8  # first nat-8 are framework, last 8 are iodine
    L=[f"&CONTROL\n  calculation = 'relax'\n  prefix = '{tag}'\n  outdir = './tmp_{tag}'",
       "  pseudo_dir = '/workspace/qe-7.5/pseudo'\n  tprnfor = .true.\n  etot_conv_thr = 1.0d-5\n  forc_conv_thr = 1.0d-3\n/",
       f"&SYSTEM\n  ibrav = 0\n  nat = {nat}\n  ntyp = {len(uniq)}\n  ecutwfc = 50\n  ecutrho = 400",
       "  occupations = 'smearing'\n  smearing = 'mv'\n  degauss = 0.01\n  nspin = 1\n  vdw_corr = 'grimme-d3'\n  dftd3_version = 4\n/",
       "&ELECTRONS\n  conv_thr = 1.0d-6\n  diagonalization = 'david'\n  diago_david_ndim = 2",
       "  mixing_mode = 'local-TF'\n  mixing_beta = 0.2\n  mixing_ndim = 24\n  electron_maxstep = 300\n/",
       "&IONS\n  ion_dynamics = 'bfgs'\n  nstep = 100\n/",
       "ATOMIC_SPECIES"]
    for e in uniq: L.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell: L.append(f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for i,(e,x,y,z) in enumerate(pos):
        flags = "0 0 0" if i < nfw else "1 1 1"   # framework frozen, iodine free
        L.append(f"  {e:<2s} {x:.10f} {y:.10f} {z:.10f}  {flags}")
    L.append("K_POINTS automatic\n  4 4 2 0 0 0\n")
    open(f"{OUT}/{tag}.in","w").write("\n".join(L)); return tag
# molecular I2 start: activated_i2_mol relax-final
cm,pm=geom_out_final(f"{ROOT}/13_Y/activated_i2_mol/loaded.out")
emit("relax442_mol_yact",cm,pm)
# polyiodide start: scf_loaded.in (Y_actc) = activated_i2 relax-final
cp,pp=geom_inp(f"{BASE}/Y_actc/scf_loaded.in")
emit("relax442_poly_yact",cp,pp)
open(f"{OUT}/relax442_manifest.txt","w").write("relax442_mol_yact\nrelax442_poly_yact\n")
print(f"wrote 2 relax inputs (framework frozen, I free, 4x4x2, 50/400)")
print(f"  molecular: nat={len(pm)}, {sum(1 for p in pm if p[0]=='I')} I atoms")
print(f"  polyiodide: nat={len(pp)}, {sum(1 for p in pp if p[0]=='I')} I atoms")
# verify if_pos flags
mol_in=open(f"{OUT}/relax442_mol_yact.in").read().splitlines()
for l in mol_in:
    if "I " in l and "1 1 1" in l: print(f"  iodine free check: {l.strip()[:60]}..."); break
for l in mol_in:
    if "C " in l and "0 0 0" in l: print(f"  framework frozen check: {l.strip()[:60]}..."); break
