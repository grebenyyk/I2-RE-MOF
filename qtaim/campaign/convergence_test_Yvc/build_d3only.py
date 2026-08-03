#!/usr/bin/env python3
"""Build PBE+D3 Gamma-point single-point inputs for every landscape config + references,
solely to extract the D3(BJ)v4 dispersion energy (geometry-only, k-independent).
Reads geometry from each <tag>.in, emits d3only_<tag>.in (PBE, Gamma, scf, grimme-d3 v4)."""
import re,os,glob
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
     "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","Tb":"Tb.pbe-n-kjpaw_psl.1.0.0.UPF",
     "Gd":"Gd.pbe-n-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"Tb":158.925,"Gd":157.250,"I":126.904}
# configs (anything with an HSE energy) + references
def has_energy(tag):
    f=f"{D}/{tag}.out"
    if not os.path.exists(f): return False
    return bool(re.search(r"^!\s*total energy\s*=",open(f,errors="ignore").read(),re.M))
tags=set()
for out in glob.glob(f"{D}/*.out"):
    t=os.path.basename(out)[:-4]
    if has_energy(t) and re.match(r'(rd_poly|rd_mol|hmol|hse_lmol|hse_lpoly)_',t): tags.add(t)
for r in ["rd_gas","rd_empty_Yst","rd_empty_Yact","rd_empty_Yvc"]: tags.add(r)
def geom(tag):
    L=open(f"{D}/{tag}.in").read().splitlines()
    ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
    ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][0]
    cell=[L[ci+1+k].split()[:3] for k in range(3)]; atoms=[];k=ai+1
    while k<len(L):
        s=L[k].split()
        if not s or s[0] in("K_POINTS","ATOMIC_SPECIES","CONSTRAINTS","OCCUPATIONS"): break
        if len(s)>=4 and re.match(r'^[A-Z]',s[0]): atoms.append((s[0],s[1],s[2],s[3]))
        k+=1
    return cell,atoms
n=0
for tag in sorted(tags):
    cell,atoms=geom(tag); uniq=[];[uniq.append(a[0]) for a in atoms if a[0] not in uniq]
    pre=f"d3only_{tag}"
    LL=[f"&CONTROL\n  calculation='scf'\n  prefix='{pre}'\n  outdir='./tmp_{pre}'",
        "  pseudo_dir='/workspace/qe-7.5/pseudo'\n  etot_conv_thr=1.0d-6\n/",
        f"&SYSTEM\n  ibrav=0\n  nat={len(atoms)}\n  ntyp={len(uniq)}\n  ecutwfc=50\n  ecutrho=400",
        "  occupations='smearing'\n  smearing='mv'\n  degauss=0.02\n  nspin=1\n  vdw_corr='grimme-d3'\n  dftd3_version=4\n/",
        "&ELECTRONS\n  conv_thr=1.0d-6\n  diagonalization='david'\n  diago_david_ndim=2\n  mixing_mode='local-TF'\n  mixing_beta=0.2\n  mixing_ndim=24\n  electron_maxstep=100\n/","ATOMIC_SPECIES"]
    for e in uniq: LL.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    LL.append("CELL_PARAMETERS angstrom")
    for r in cell: LL.append(f"  {r[0]} {r[1]} {r[2]}")
    LL.append("ATOMIC_POSITIONS crystal")
    for e,x,y,z in atoms: LL.append(f"  {e:<2s} {x} {y} {z}")
    LL.append("K_POINTS automatic\n  1 1 1 0 0 0\n")
    open(f"{D}/{pre}.in","w").write("\n".join(LL)); n+=1
print(f"built {n} d3only inputs (configs + references)")
