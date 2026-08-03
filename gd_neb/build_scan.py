#!/usr/bin/env python3
"""Constrained-coordinate scan: ch1 transformation (2 I2 -> I4 unit) by constraining
the I1-I4 distance (primary chain-junction former) from molecular (4.17 A) to chain
(3.14 A). Framework + ch2 iodines FROZEN (reference state); ch1 I atoms (I1,I3,I4,I5)
relaxable; I1-I4 distance constrained to target d. Each window = a constrained MIN =>
clean SCF, no NEB force-feedback disease. Usage: build_scan.py <d_angstrom> <label>"""
import sys, os, numpy as np, ase.io
d_A = float(sys.argv[1]); label = sys.argv[2]
d_bohr = d_A / 0.5291772107   # QE constr_target is in bohr (a.u.)
BASE = "/workspace/dimitry/calc/qe/i2_mof"; GD = f"{BASE}/gd_neb/scan_{label}"
mol = ase.io.read(f"{BASE}/08_Gd/activated_i2/POSCAR_final_optimized")
cell = mol.cell.array
sym = np.array(mol.get_chemical_symbols())
# ch1 I atoms (mol numbering): I1,I3,I4,I5 -> relaxable. Everything else frozen.
ch1_I = {"I": [0, 2, 3, 4]}  # 0-based indices among I atoms
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF",
     "O":"O.pbe-kjpaw.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"Gd":157.25,"H":1.008,"O":15.999,"I":126.904}
uniq=[]; [uniq.append(s) for s in sym if s not in uniq]
I_global = [i+1 for i,s in enumerate(sym) if s=="I"]  # 1-based global indices of I atoms
I1g = I_global[0]; I4g = I_global[3]   # I1 and I4 global atom numbers
os.makedirs(GD, exist_ok=True)
L = ["&CONTROL", f"  calculation = 'relax'", f"  prefix = 'scan_{label}'",
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
ii = -1
for s,p in zip(sym, mol.positions):
    if s == "I":
        ii += 1
        relaxable = ii in ch1_I["I"]
        ifp = "1 1 1" if relaxable else "0 0 0"
    else:
        ifp = "0 0 0"
    L.append(f"{s:<3s} {p[0]: .10f} {p[1]: .10f} {p[2]: .10f} {ifp}")
# CONSTRAINTS card: 1 constraint, distance I1-I4 = target (bohr)
L += ["CONSTRAINTS", f"  1  1.0d-8", f"  'distance'  {I1g}  {I4g}  {d_bohr:.6f}",
      "K_POINTS gamma", ""]
open(f"{GD}/relax.in","w").write("\n".join(L))
print(f"scan_{label}: constrain I1(atom{I1g})-I4(atom{I4g}) = {d_A:.2f} A ({d_bohr:.4f} bohr); ch1 relax, rest frozen -> {GD}/relax.in")
