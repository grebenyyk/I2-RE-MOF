#!/usr/bin/env python3
"""Build HSE scf + projwfc inputs for the 3 Y-MOF variants (st/act/vc) to obtain
hybrid-level PDOS and confirm the I->framework charge transfer is not a PBE
self-interaction artifact. Geometry is read from each variant's PBE scf_loaded.in
(the SAME structures used for the PBE pdos_Y_3variants plot -> apples-to-apples).

Recipe mirrors (a) the converged Y_str/hse/scf.in (HSE, 50/400, mv/0.03, ACE,
33 iters ~10 min) and (b) the WORKING PBE proj.in (ngauss=0, degauss=0.01,
DeltaE=0.01). The only addition is wf_collect=.true. so projwfc can read wfc1.dat.

Previous attempt failed because the scf was KILLED mid post-convergence EXX pass
(before wfc1.dat was written); running to JOB DONE with wf_collect=.true. fixes it.
"""
import os, re
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
SYS=[("Y_str","st"),("Y_actc","act"),("Y_deep","vc")]
PPMAP={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
       "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"I":126.904}

def parse_loaded(path):
    t=open(path,errors="ignore").read()
    m=re.search(r"CELL_PARAMETERS[^\n]*\n((?:[^\n]+\n){3})",t)
    n=re.findall(r"-?\d+\.?\d*",m.group(1))
    cell=[[float(n[i*3+j]) for j in range(3)] for i in range(3)]
    ms=re.search(r"ATOMIC_SPECIES\n((?:[A-Z][a-z]?.*\n)+)",t)
    spec=[ln.split()[0] for ln in ms.group(1).split("\n") if ln.strip() and ln.split()[0] in PPMAP]
    mp=re.search(r"ATOMIC_POSITIONS[^\n]*\n((?:.+\n)+?)(?:K_POINTS|\Z)",t)
    atoms=[]
    for ln in mp.group(1).split("\n"):
        mm=re.match(r"\s*([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",ln)
        if mm and mm.group(1) in MASS: atoms.append((mm.group(1),[float(mm.group(i)) for i in (2,3,4)]))
    return cell,spec,atoms

def write_scf(path,prefix,cell,spec,atoms):
    L=["&CONTROL","  calculation = 'scf'",f"  prefix = '{prefix}'","  outdir = './tmp_hse'",
       "  pseudo_dir = '/workspace/qe-7.5/pseudo'","  verbosity = 'high'",
       "  wf_collect = .true.","  disk_io = 'high'","/",
       "&SYSTEM","  ibrav = 0",f"  nat = {len(atoms)}",f"  ntyp = {len(spec)}",
       "  ecutwfc = 50","  ecutrho = 400","  occupations = 'smearing'","  smearing = 'mv'",
       "  degauss = 0.03","  nspin = 1","  input_dft = 'hse'",
       "  exxdiv_treatment = 'gygi-baldereschi'","  x_gamma_extrapolation = .true.",
       "  exx_fraction = 0.25","  screening_parameter = 0.106",
       "  ecutfock = 150","/",
       "&ELECTRONS","  conv_thr = 1.0d-6","  diagonalization = 'david'","  diago_david_ndim = 2",
       "  mixing_mode = 'local-TF'","  mixing_beta = 0.1","  mixing_ndim = 20",
       "  electron_maxstep = 600","/","ATOMIC_SPECIES"]
    for s in spec: L.append(f"{s:<3s} {MASS[s]:>10.3f}  {PPMAP[s]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell: L.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for s,fr in atoms: L.append(f"{s:<3s} {fr[0]: .12f} {fr[1]: .12f} {fr[2]: .12f}")
    L+=["K_POINTS gamma",""]
    open(path,"w").write("\n".join(L))

def write_proj(path,prefix):
    open(path,"w").write(
        f"&PROJWFC\n  prefix = '{prefix}'\n  outdir = './tmp_hse'\n  filpdos = 'hse_pdos'\n"
        "  DeltaE = 0.01\n  ngauss = 0\n  degauss = 0.01\n/\n")

for d,tag in SYS:
    cell,spec,atoms=parse_loaded(f"{BASE}/{d}/scf_loaded.in")
    assert len(atoms)==80, f"{d}: {len(atoms)} atoms (expected 80)"
    assert set(spec)=={"C","H","O","Y","I"}, f"{d}: species {spec}"
    os.makedirs(f"{BASE}/{d}/hse",exist_ok=True)
    prefix=f"{d}_hse"
    write_scf(f"{BASE}/{d}/hse/scf.in",prefix,cell,spec,atoms)
    write_proj(f"{BASE}/{d}/hse/proj.in",prefix)
    print(f"{d} ({tag}): {len(atoms)} atoms, {len(spec)} species -> hse/scf.in + hse/proj.in  prefix={prefix}")
print("done")
