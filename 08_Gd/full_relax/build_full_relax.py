#!/usr/bin/env python3
"""
Build QE inputs for a FULL-framework relaxation of Gd I2@MOF (ALL atoms movable,
fixed cell). Tests whether the frozen-framework E_int/I2 = -2.0753 eV is an artifact
of freezing the framework.

Writes:
  empty_full.in  : bare Gd framework (72 atoms, I removed), ALL atoms movable, relax
  loaded_full.in : full loaded structure (80 atoms from uniform_chain/POSCAR_final),
                   ALL atoms movable, relax

NO if_pos freezing anywhere (every atom moves). Compare with build_inputs.py, which
sets if_pos 0 0 0 for framework atoms -- that is deliberately NOT done here.

E_int formula once results exist:
  E_int/I2 = (E_loaded_full - E_empty_full - 4*(-757.80004425)) * 13.605693009 / 4  [eV]
           then * 96.485 for kJ/mol
"""
import os
import numpy as np
import ase.io

PSEUDO_DIR = "/workspace/qe-7.5/pseudo"
ECUTWFC = 50
ECUTRHO = 400
DEGAUSS = "0.01"
MIXING_BETA = "0.2"

# pseudos + masses
SPECIES = {
    "C":  ("C.pbe-n-kjpaw_psl.1.0.0.UPF", 12.011),
    "H":  ("H.pbe-kjpaw.UPF",              1.008),
    "O":  ("O.pbe-kjpaw.UPF",             15.999),
    "I":  ("I.pbe-n-kjpaw_psl.1.0.0.UPF", 126.904),
    "Gd": ("Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF", 157.25),
}

POSCAR = "/workspace/dimitry/calc/qe/i2_mof/08_Gd/uniform_chain/POSCAR_final"
OUTDIR = "/workspace/dimitry/calc/qe/i2_mof/08_Gd/full_relax"


def species_table(symbols):
    uniq = []
    for s in symbols:
        if s not in uniq:
            uniq.append(s)
    return uniq


def write_input(path, syms_list, frac_list, prefix):
    cell = np.array(ase.io.read(POSCAR).cell)
    uniq = species_table(syms_list)
    nat = len(syms_list)
    ntyp = len(uniq)
    lines = []
    lines += [
        "&CONTROL",
        "  calculation = 'relax'",
        f"  prefix = '{prefix}'",
        "  outdir = './tmp'",
        f"  pseudo_dir = '{PSEUDO_DIR}'",
        "  tprnfor = .true.",
        "  tstress = .true.",
        "  etot_conv_thr = 1.0d-4",
        "  forc_conv_thr = 1.0d-2",
        "/",
        "&SYSTEM",
        "  ibrav = 0",
        f"  nat = {nat}",
        f"  ntyp = {ntyp}",
        f"  ecutwfc = {ECUTWFC}",
        f"  ecutrho = {ECUTRHO}",
        "  occupations = 'smearing'",
        "  smearing = 'mv'",
        f"  degauss = {DEGAUSS}",
        "  nspin = 1",
        "  vdw_corr = 'grimme-d3'",
        "  dftd3_version = 4",
        "/",
        "&ELECTRONS",
        "  conv_thr = 1.0d-6",
        "  diagonalization = 'david'",
        "  diago_david_ndim = 4",
        "  mixing_mode = 'local-TF'",
        f"  mixing_beta = {MIXING_BETA}",
        "  mixing_ndim = 16",
        "  electron_maxstep = 300",
        "/",
        "&IONS",
        "  ion_dynamics = 'bfgs'",
        "/",
        "ATOMIC_SPECIES",
    ]
    for el in uniq:
        p, m = SPECIES[el]
        lines.append(f"{el:<3s} {m:>10.3f}  {p}")
    lines += ["CELL_PARAMETERS angstrom"]
    for r in cell:
        lines.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}")
    # ALL atoms movable: no if_pos column at all
    lines += ["ATOMIC_POSITIONS crystal"]
    for el, f in zip(syms_list, frac_list):
        lines.append(f"{el:<3s} {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}")
    lines += ["K_POINTS automatic", "  1 1 1 0 0 0", ""]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
    return path


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    atoms = ase.io.read(POSCAR)
    frac = atoms.get_scaled_positions(wrap=False)
    syms = list(atoms.get_chemical_symbols())
    assert len(syms) == 80, f"expected 80 atoms in loaded POSCAR, got {len(syms)}"

    # empty: drop I atoms, keep order
    fw_idx = [i for i, s in enumerate(syms) if s != "I"]
    fw_syms = [syms[i] for i in fw_idx]
    fw_frac = frac[fw_idx]
    assert len(fw_syms) == 72, f"expected 72 framework atoms, got {len(fw_syms)}"

    # loaded: all 80 atoms
    ld_syms = syms
    ld_frac = frac

    ep = write_input(os.path.join(OUTDIR, "empty_full.in"), fw_syms, fw_frac, "Gd_empty_full")
    lp = write_input(os.path.join(OUTDIR, "loaded_full.in"), ld_syms, ld_frac, "Gd_loaded_full")
    print(f"empty  ({len(fw_syms)} atoms) -> {ep}")
    print(f"loaded ({len(ld_syms)} atoms) -> {lp}")


if __name__ == "__main__":
    main()
