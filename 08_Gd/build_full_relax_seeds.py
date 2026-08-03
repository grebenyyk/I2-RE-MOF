#!/usr/bin/env python3
"""
Build TWO new "full-relax" Gd I2@MOF sorption inputs (cell FIXED, ALL atoms movable).

Seeds:
  full_relax_chain/relax.in  <- uniform_chain/POSCAR_final       (chain, ~ -200 kJ/mol)
  full_relax_mol/relax.in    <- activated_i2/POSCAR_final_optimized (molecular I2, ~ -64.5 kJ/mol)

Namelist mirrors 08_Gd/uniform_chain/loaded.in (ecutwfc=50, ecutrho=400, PBE+D3(BJ),
Gamma-only K_POINTS, bfgs) with the FULL-RELAX changes:
  - calculation = 'relax'
  - every atom movable (NO if_pos column anywhere -- framework + iodine both free)
  - cell fixed (no vc-relax)
  - Gd-SCF-safe recipe: degauss=0.01 (mv), mixing_beta=0.2 (PROVEN on this exact
    chain seed -- uniform_chain/loaded.out converged first SCF in 70 iters with
    beta=0.2; the earlier attempt here with beta=0.1 hit the 300-iter SCF cap and
    failed), diago_david_ndim=2 (VRAM-safe), electron_maxstep=300, nstep=200,
    conv_thr=1e-5 (matches the seed loaded.in)
  - UNIQUE prefix/outdir per calc (Gd_fr_chain / Gd_fr_mol) to avoid .save collisions

Does NOT touch any existing campaign file. Does NOT overwrite full_relax/.
"""
import os
import numpy as np
import ase.io

PSEUDO_DIR = "/workspace/qe-7.5/pseudo"
BASE = "/workspace/dimitry/calc/qe/i2_mof/08_Gd"

# pseudos + masses (EXACT names -- do NOT glob)
SPECIES = {
    "C":  ("C.pbe-n-kjpaw_psl.1.0.0.UPF",  12.011),
    "H":  ("H.pbe-kjpaw.UPF",                1.008),
    "O":  ("O.pbe-kjpaw.UPF",               15.999),
    "I":  ("I.pbe-n-kjpaw_psl.1.0.0.UPF",  126.904),
    "Gd": ("Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF", 157.250),
}

# Seed POSCARs (fractional/direct coords). Both have the SAME unit cell as the
# 08_Gd campaign (9.0752 / 9.6427 / c-vector with -3.898... 16.784...).
SEEDS = [
    ("full_relax_chain", "Gd_fr_chain",
     os.path.join(BASE, "uniform_chain/POSCAR_final")),
    ("full_relax_mol",   "Gd_fr_mol",
     os.path.join(BASE, "activated_i2/POSCAR_final_optimized")),
]


def species_order(symbols):
    """Order matching the campaign ATOMIC_SPECIES block: C Gd H O I."""
    preferred = ["C", "Gd", "H", "O", "I"]
    seen = []
    for s in symbols:
        if s not in seen:
            seen.append(s)
    # sort by preferred order if all present there, else keep first-seen
    if all(s in preferred for s in seen):
        return [s for s in preferred if s in seen]
    return seen


def write_input(atoms, out_dir, prefix):
    os.makedirs(out_dir, exist_ok=True)
    cell = np.array(atoms.get_cell())
    syms = list(atoms.get_chemical_symbols())
    frac = atoms.get_scaled_positions(wrap=False)
    nat = len(syms)
    uniq = species_order(syms)
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
        "  etot_conv_thr = 1.0d-5",
        "  forc_conv_thr = 1.0d-3",
        "  nstep = 200",
        "/",
        "&SYSTEM",
        "  ibrav = 0",
        f"  nat = {nat}",
        f"  ntyp = {ntyp}",
        "  ecutwfc = 50",
        "  ecutrho = 400",
        "  occupations = 'smearing'",
        "  smearing = 'mv'",
        "  degauss = 0.01",
        "  nspin = 1",
        "  vdw_corr = 'grimme-d3'",
        "  dftd3_version = 4",
        "/",
        "&ELECTRONS",
        "  conv_thr = 1.0d-5",
        "  diagonalization = 'david'",
        "  diago_david_ndim = 2",
        "  mixing_mode = 'local-TF'",
        "  mixing_beta = 0.2",
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
    # ALL atoms movable: NO if_pos column anywhere (framework + I both free).
    lines += ["ATOMIC_POSITIONS crystal"]
    for el, f in zip(syms, frac):
        lines.append(f"{el:<3s} {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}")
    lines += ["K_POINTS automatic", "  1 1 1 0 0 0", ""]

    inp_path = os.path.join(out_dir, "relax.in")
    with open(inp_path, "w") as fh:
        fh.write("\n".join(lines))

    # sanity: verify no atom has the same fractional coords as another within 1e-9
    frac_arr = np.array(frac)
    for i in range(nat):
        d = np.abs(frac_arr - frac_arr[i])
        d = np.minimum(d, 1.0 - d)
        for j in range(i + 1, nat):
            if np.all(d[j] < 1e-9):
                raise RuntimeError(f"atoms {i}({syms[i]}) and {j}({syms[j]}) coincide")
    return inp_path, nat


def main():
    for subdir, prefix, poscar in SEEDS:
        if not os.path.isfile(poscar):
            raise FileNotFoundError(poscar)
        atoms = ase.io.read(poscar)
        out_dir = os.path.join(BASE, subdir)
        inp_path, nat = write_input(atoms, out_dir, prefix)
        cell = atoms.get_cell()
        print(f"{subdir}: {nat} atoms, prefix={prefix}")
        print(f"  seed   : {poscar}")
        print(f"  input  : {inp_path}")
        print(f"  cell a: {cell[0]}")
        print(f"  cell b: {cell[1]}")
        print(f"  cell c: {cell[2]}")
        # quick atom-count breakdown
        from collections import Counter
        cnt = Counter(atoms.get_chemical_symbols())
        print(f"  counts : {dict(cnt)}")


if __name__ == "__main__":
    main()
