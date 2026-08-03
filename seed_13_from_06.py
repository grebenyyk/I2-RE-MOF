#!/usr/bin/env python3
"""Seed 13_Y with 06_Tb's final relaxed I2 positions (isostructural pts frameworks,
related by origin shift Y = Tb + (0, 0.5, 0.5) mod 1), to help 13_Y escape its plateau.

Reads 06_Tb's final I2 fractional coords from 06_Tb/POSCAR_final (regenerate first via
`python qe_out_to_poscar.py 06_Tb/loaded_tight.out`), shifts them, and writes
13_Y/loaded_seeded.in: 13_Y framework (frozen) + shifted I2 (movable), tight settings.
"""
import numpy as np, ase.io

PSEUDO_DIR = "/workspace/qe-7.5/pseudo"
SHIFT = np.array([0.0, 0.5, 0.5])  # Y = Tb + shift (mod 1)
PSEUDOS = {"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
           "N":"N.pbe-n-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF",
           "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF"}
MASS = {"C":12.011,"H":1.008,"O":15.999,"N":14.007,"I":126.904,"Y":88.906}

# 13_Y framework (frozen) from the original POSCAR
Y = ase.io.read("poscars_unitcell/POSCAR_13_2020_Y_pts_3_ASR_1_I2loaded_1000Pa_UC")
ycell = np.array(Y.cell)
yf = Y.get_scaled_positions(wrap=False)
ysym = list(Y.get_chemical_symbols())
fw_idx = [i for i,s in enumerate(ysym) if s != "I"]

# 06_Tb final relaxed I2 (from its POSCAR_final)
Tb = ase.io.read("06_Tb/POSCAR_final")
tf = Tb.get_scaled_positions(wrap=False)
tsy = Tb.get_chemical_symbols()
tb_I = [(tf[i] + SHIFT) % 1.0 for i,s in enumerate(tsy) if s == "I"]
print(f"06_Tb final I atoms: {len(tb_I)}  -> shifted into 13_Y frame")

# assemble: framework (frozen) + shifted I2 (movable)
syms = [ysym[i] for i in fw_idx] + ["I"]*len(tb_I)
frac = [yf[i] for i in fw_idx] + tb_I
cell = ycell
uniq = []
for s in syms:
    if s not in uniq: uniq.append(s)

L = ["&CONTROL","  calculation = 'relax'","  prefix = '13_Y_seeded'","  outdir = './tmp'",
     f"  pseudo_dir = '{PSEUDO_DIR}'","  tprnfor = .true.","  tstress = .true.",
     "  nstep = 200","  etot_conv_thr = 1.0d-5","  forc_conv_thr = 1.0d-3","/",
     "&SYSTEM","  ibrav = 0",f"  nat = {len(syms)}",f"  ntyp = {len(uniq)}",
     "  ecutwfc = 50","  ecutrho = 400","  occupations = 'smearing'","  smearing = 'mv'",
     "  degauss = 0.01","  nspin = 1","  vdw_corr = 'grimme-d3'","  dftd3_version = 4","/",
     "&ELECTRONS","  conv_thr = 1.0d-6","  diagonalization = 'david'","  diago_david_ndim = 4",
     "  mixing_mode = 'local-TF'","  mixing_beta = 0.2","  mixing_ndim = 16","  electron_maxstep = 300","/",
     "&IONS","  ion_dynamics = 'bfgs'","/","ATOMIC_SPECIES"]
for s in uniq: L.append(f"{s:<3s} {MASS[s]:>10.3f}  {PSEUDOS[s]}")
L.append("CELL_PARAMETERS angstrom")
for r in cell: L.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}")
L.append("ATOMIC_POSITIONS crystal")
for s,f in zip(syms,frac):
    L.append(f"{s:<3s} {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}  " + ("1 1 1" if s=="I" else "0 0 0"))
L += ["K_POINTS automatic","  1 1 1 0 0 0",""]
open("13_Y/loaded_seeded.in","w").write("\n".join(L))
print(f"wrote 13_Y/loaded_seeded.in: {len(syms)} atoms ({len(fw_idx)} fw + {len(tb_I)} I)")
