#!/usr/bin/env python3
"""Collective constrained scan: drive the GLOBAL I-I coordination number (type_coord)
from molecular (~1.0) to chain (~2.0). This forces ALL I-I contacts to form collectively
across both channels, so the system cannot route around the barrier by relaxing one
channel (the flaw of the single-distance scan). Framework frozen, all 8 I relaxable.
Usage: build_cscan.py <target_CN> <label>"""
import sys, os, numpy as np, ase.io
cn = float(sys.argv[1]); label = sys.argv[2]
RC_BOHR = 6.05      # 3.2 A cutoff (between I2 bond 2.7 and chain contact ~3.1)
SMOOTH = 0.2        # Fermi-Dirac smoothing (bohr)
BASE = "/workspace/dimitry/calc/qe/i2_mof"; GD = f"{BASE}/gd_neb/cscan_{label}"
mol = ase.io.read(f"{BASE}/08_Gd/activated_i2/POSCAR_final_optimized")
cell = mol.cell.array; sym = np.array(mol.get_chemical_symbols())
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF",
     "O":"O.pbe-kjpaw.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"Gd":157.25,"H":1.008,"O":15.999,"I":126.904}
uniq=[]; [uniq.append(s) for s in sym if s not in uniq]
I_type = uniq.index("I") + 1   # 1-based type index for I
os.makedirs(GD, exist_ok=True)
L = ["&CONTROL", "  calculation = 'relax'", f"  prefix = 'cscan_{label}'",
     f"  outdir = '{GD}/tmp'", "  pseudo_dir = '/workspace/qe-7.5/pseudo'", "  tprnfor = .true.",
     "  nstep = 100","  dt = 20.0","  etot_conv_thr = 1.0d-5","/",
     "&SYSTEM","  ibrav = 0",f"  nat = {len(sym)}",f"  ntyp = {len(uniq)}","  ecutwfc = 50","  ecutrho = 400",
     "  occupations = 'smearing'","  smearing = 'mv'","  degauss = 0.01","  nspin = 1",
     "  vdw_corr = 'grimme-d3'","  dftd3_version = 4","/",
     "&ELECTRONS","  conv_thr = 1.0d-4","  diagonalization = 'david'","  diago_david_ndim = 4",
     "  mixing_mode = 'local-TF'","  mixing_beta = 0.1","  mixing_ndim = 16","  electron_maxstep = 200","/",
     "&IONS","  ion_dynamics = 'damp'","/",
     "ATOMIC_SPECIES"]
for s in uniq: L.append(f"{s:<3s} {MASS[s]:>10.3f}  {PSE[s]}")
L += ["CELL_PARAMETERS angstrom"]
for r in cell: L.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}")
L += ["ATOMIC_POSITIONS angstrom"]
for s,p in zip(sym, mol.positions):
    ifp = "0 0 0" if s != "I" else "1 1 1"   # framework frozen, ALL iodines relax
    L.append(f"{s:<3s} {p[0]: .10f} {p[1]: .10f} {p[2]: .10f} {ifp}")
# type_coord constraint: global I-I coordination = target CN
L += ["CONSTRAINTS", "  1  1.0d-8",
      f"  'type_coord'  {I_type}  {I_type}  {RC_BOHR}  {SMOOTH}  {cn}",
      "K_POINTS gamma", ""]
open(f"{GD}/relax.in","w").write("\n".join(L))
print(f"cscan_{label}: constrain I-I type_coord(I_type={I_type}, rc=3.2A, smooth=0.2) = {cn}; fw frozen, all I relax -> {GD}/relax.in")
