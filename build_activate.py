#!/usr/bin/env python3
"""
Stage-1 ACTIVATE: build a bare-framework (I removed) RELAX input for a MOF system.

Mirrors 08_Gd/activate_empty/relax.in exactly:
  - calculation = 'relax', ALL atoms movable (NO if_pos), fixed cell (ibrav=0)
  - nstep = 200  (framework relaxation needs many steps; default 50 caps too early)
  - conv_thr = 1e-6, etot_conv_thr = 1e-4, forc_conv_thr = 1e-2
  - degauss = 0.01, smearing 'mv'  (bare frameworks are gapped)
  - PBE + D3(BJ) dftd3_version=4, 50/400 PAW, Gamma point, nspin=1, no Hubbard
  - davidson / local-TF mixing beta 0.2 ndim 16, electron_maxstep 300
  - ions bfgs

Atom order = the POSCAR order with I atoms dropped (matches the Gd template:
C Gd H O, i.e. POSCAR species order minus I).  Prefix = 'activate_empty'.

Usage:
  python build_activate.py <system>            # uses poscars_unitcell/POSCAR_<...>
  python build_activate.py <system> <POSCAR>   # explicit POSCAR
  python build_activate.py --all               # build the 7 non-Gd systems
  python build_activate.py --verify 08_Gd      # check existing matches template

After this runs, pw.x produces 08_Gd/activate_empty/relax.out with the final
activated-framework ATOMIC_POSITIONS, which build_activated_i2.py reads.
"""
import os, sys
import numpy as np
import ase.io

PSEUDO_DIR = "/workspace/qe-7.5/pseudo"
ECUTWFC, ECUTRHO = 50, 400
DEGAUSS = 0.01

# (pseudo, mass) -- PBE PAW PSlibrary 1.0.0, identical to build_inputs.py
COMMON = {
    "C": ("C.pbe-n-kjpaw_psl.1.0.0.UPF", 12.011),
    "H": ("H.pbe-kjpaw.UPF", 1.008),
    "O": ("O.pbe-kjpaw.UPF", 15.999),
    "N": ("N.pbe-n-kjpaw_psl.1.0.0.UPF", 14.007),
    "I": ("I.pbe-n-kjpaw_psl.1.0.0.UPF", 126.904),
}
METALS = {
    "Y":  ("Y.pbe-spn-kjpaw_psl.1.0.0.UPF",  88.906),
    "Cd": ("Cd.pbe-n-kjpaw_psl.1.0.0.UPF",  112.414),
    "Gd": ("Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF", 157.25),
    "Tb": ("Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF", 158.925),
}

# system -> POSCAR filename in poscars_unitcell/  (same map as build_inputs.py)
POSCARS = {
    "13_Y":  "POSCAR_13_2020_Y_pts_3_ASR_1_I2loaded_1000Pa_UC",
    "06_Tb": "POSCAR_06_2019_Tb_pts_3_ASR_1_I2loaded_1000Pa_UC",
    "08_Gd": "POSCAR_08_2020_Gd_pts_3_ASR_1_I2loaded_1000Pa_UC",
    "12_Tb": "POSCAR_12_2020_Tb_nan_3_ASR_7_I2loaded_1000Pa_UC",
    "05_Cd": "POSCAR_05_oc8b00378_si_004_ASR_pacman_I2loaded_1000Pa_UC",
    "07_Cd": "POSCAR_07_oc8b00378_si_006_ASR_pacman_I2loaded_1000Pa_UC",
    "10_Cd": "POSCAR_10_oc8b00378_si_003_ASR_pacman_I2loaded_1000Pa_UC",
    "17_Cd": "POSCAR_17_2016_Cd_nan_3_ASR_1_I2loaded_1000Pa_UC",
}


def detect_metal(symbols):
    metals = [s for s in symbols if s not in COMMON]
    assert len(set(metals)) == 1, f"expected one metal, got {set(metals)}"
    return metals[0]


def pseudo_mass(el, metal):
    if el == metal:
        return METALS[metal]
    return COMMON[el]


def species_table(symbols):
    """Ordered unique species list (order of first appearance)."""
    uniq = []
    for s in symbols:
        if s not in uniq:
            uniq.append(s)
    return uniq


def build(poscar, outdir, label):
    """Write <outdir>/relax.in : bare-framework (I removed) all-movable relax."""
    atoms = ase.io.read(poscar)
    cell = np.array(atoms.cell)
    frac = atoms.get_scaled_positions(wrap=False)
    syms = list(atoms.get_chemical_symbols())
    metal = detect_metal(syms)

    # drop I atoms, keep POSCAR order
    keep = [i for i, s in enumerate(syms) if s != "I"]
    fw_syms = [syms[i] for i in keep]
    fw_frac = frac[keep]
    nat = len(fw_syms)

    uniq = species_table(fw_syms)
    ntyp = len(uniq)

    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "relax.in")

    L = []
    L += [
        "&CONTROL",
        "  calculation = 'relax'",
        "  prefix = 'activate_empty'",
        "  outdir = './tmp'",
        f"  pseudo_dir = '{PSEUDO_DIR}'",
        "  tprnfor = .true.",
        "  tstress = .true.",
        "  etot_conv_thr = 1.0d-4",
        "  forc_conv_thr = 1.0d-2",
        "  nstep = 200",
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
        p, m = pseudo_mass(el, metal)
        L.append(f"{el:<3s} {m:>10.3f}  {p}")
    L += ["CELL_PARAMETERS angstrom"]
    for r in cell:
        L.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}")
    L += ["ATOMIC_POSITIONS crystal"]
    # ALL atoms movable -> no if_pos (matches the Gd template exactly)
    for el, f in zip(fw_syms, fw_frac):
        L.append(f"{el:<3s} {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}")
    L += ["K_POINTS automatic", "  1 1 1 0 0 0", ""]

    with open(path, "w") as fh:
        fh.write("\n".join(L))
    print(f"[{label}] metal={metal} nat={nat} (loaded {nat+syms.count('I')} "
          f"minus {syms.count('I')} I) -> {path}")
    return path


def main():
    base = "/workspace/dimitry/calc/qe/i2_mof"
    posdir = os.path.join(base, "poscars_unitcell")
    args = sys.argv[1:]
    if not args:
        print(__doc__); sys.exit(1)
    if args[0] == "--all":
        # 08_Gd already exists -- skip, do not overwrite
        for sys_ in ["13_Y", "06_Tb", "12_Tb", "05_Cd", "07_Cd", "10_Cd", "17_Cd"]:
            build(os.path.join(posdir, POSCARS[sys_]),
                  os.path.join(base, sys_, "activate_empty"), sys_)
        return
    if args[0] == "--verify":
        # confirm the existing 08_Gd/activate_empty/relax.in matches template settings
        sys_ = args[1]
        p = os.path.join(base, sys_, "activate_empty", "relax.in")
        with open(p) as fh:
            t = fh.read()
        checks = {
            "calculation = 'relax'": "calc = relax",
            "nstep = 200": "nstep 200",
            "etot_conv_thr = 1.0d-4": "etot 1e-4",
            "forc_conv_thr = 1.0d-2": "forc 1e-2",
            "conv_thr = 1.0d-6": "conv 1e-6",
            "degauss = 0.01": "degauss 0.01",
            "nspin = 1": "nspin 1",
            "vdw_corr = 'grimme-d3'": "vdw",
            "dftd3_version = 4": "d3v4",
            "diagonalization = 'david'": "david",
            "mixing_mode = 'local-TF'": "local-TF",
            "mixing_beta = 0.2": "beta 0.2",
            "ecutwfc = 50": "ecutwfc 50",
            "ecutrho = 400": "ecutrho 400",
        }
        bad = [name for s, name in checks.items() if s not in t]
        print(f"[{sys_}] activate_empty/relax.in: ", end="")
        if bad:
            print(f"MISSING -> {bad}")
        else:
            print("all template settings present OK")
        return
    sys_ = args[0]
    poscar = args[1] if len(args) > 1 else os.path.join(posdir, POSCARS[sys_])
    build(poscar, os.path.join(base, sys_, "activate_empty"), sys_)


if __name__ == "__main__":
    main()
