#!/usr/bin/env python3
"""
Stage-2 I2@ACTIVATED: place I2 (seeded from the strained-cell final geometry)
on the ACTIVATED framework and relax ONLY the I atoms.

Mirrors 08_Gd/activated_i2/loaded.in exactly:
  - calculation = 'relax'
  - framework FROZEN (if_pos 0 0 0) at the activated geometry
  - I atoms MOVABLE (if_pos 1 1 1), seeded from the strained-cell final I2 coords
  - conv_thr = 1e-5, etot_conv_thr = 1e-5, forc_conv_thr = 1e-3
  - nstep default (50) -- only the light I atoms move, converges fast
  - degauss per-system: 0.01 (rare-earth pts/nan), 0.005 (small-gap Cd)
  - PBE + D3(BJ) dftd3_version=4, 50/400 PAW, Gamma, nspin=1, no Hubbard
  - davidson / local-TF mixing beta 0.2 ndim 16, electron_maxstep 300
  - prefix = '<system>_activated_i2_loaded'

Inputs this script needs:
  (1) activated framework geometry from <system>/activate_empty/relax.out
      (the LAST "ATOMIC_POSITIONS (crystal)" block, i.e. the converged relaxed
       framework). This file only exists AFTER stage-1 (build_activate.py +
       pw.x) has run.
  (2) I2 seed fractional coords from i2_seed_coords.json (built once from the
      strained-cell POSCAR_final files; cell is fixed so fractional transfers).

The cell is the SAME fixed cell as activation (fixed-cell throughout), so the I
fractional coords transfer directly from the strained-cell campaign. Framework
atoms are reordered to match the relax.out species order; I atoms are appended
last (movable).

Usage:
  python build_activated_i2.py <system>            # build one
  python build_activated_i2.py --all               # build all 8 (needs all relax.out)
  python build_activated_i2.py --check              # report which relax.out exist

Interaction energy (apply once both stages + the loaded relax finish):
  E_int/I2 = (E_loaded_activated - E_empty_activated - N_I2*(-757.80004425))
             * 13.605693009 / N_I2   [eV]
           (* 96.485 for kJ/mol)
where E_empty_activated = final total energy of activate_empty/relax.out and
      E_loaded_activated = final total energy of activated_i2/loaded.out.
E(I2_vac) = -757.80004425 Ry (shared reference, PBE+D3(BJ) 50/400).
"""
import os, sys, re, json
import numpy as np

PSEUDO_DIR = "/workspace/qe-7.5/pseudo"
ECUTWFC, ECUTRHO = 50, 400
E_I2_VAC = -757.80004425

# per-system stage-2 degauss (matches the strained-cell campaign)
DEGAUSS = {
    "13_Y": 0.01, "06_Tb": 0.01, "08_Gd": 0.01, "12_Tb": 0.01,
    "05_Cd": 0.005, "07_Cd": 0.005, "10_Cd": 0.005, "17_Cd": 0.005,
}
SYSTEMS = list(DEGAUSS.keys())

# (pseudo, mass) -- PBE PAW PSlibrary 1.0.0
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


def detect_metal(symbols):
    metals = [s for s in symbols if s not in COMMON]
    assert len(set(metals)) == 1, f"expected one metal, got {set(metals)}"
    return metals[0]


def pseudo_mass(el, metal):
    return METALS[metal] if el == metal else COMMON[el]


def species_table(symbols):
    uniq = []
    for s in symbols:
        if s not in uniq:
            uniq.append(s)
    return uniq


def read_cell_from_input(relax_in):
    """Read the CELL_PARAMETERS (angstrom) block from a pw.x input. The cell is
    FIXED during activation (ibrav=0, no cell relaxation), so this is identical
    to the final cell in relax.out (pw.x does not echo a CELL_PARAMETERS card
    for a fixed-cell run)."""
    with open(relax_in) as f:
        lines = f.readlines()
    cp_idx = [i for i, l in enumerate(lines)
              if "CELL_PARAMETERS" in l and "angstrom" in l]
    if not cp_idx:
        raise RuntimeError(f"no CELL_PARAMETERS (angstrom) in {relax_in}")
    i = cp_idx[-1]
    return np.array([[float(x) for x in lines[i + 1 + k].split()[:3]]
                     for k in range(3)])


def read_activated_framework(relax_out):
    """Parse the LAST 'ATOMIC_POSITIONS (crystal)' block from a pw.x relax output.

    Returns (symbols, frac[N,3]) -- the converged activated-framework geometry.
    The cell is read separately from the input (fixed-cell run; see
    read_cell_from_input).
    """
    with open(relax_out) as f:
        lines = f.readlines()
    ap_idx = [i for i, l in enumerate(lines) if "ATOMIC_POSITIONS" in l]
    if not ap_idx:
        raise RuntimeError(f"no ATOMIC_POSITIONS in {relax_out}")
    i = ap_idx[-1]
    syms, frac = [], []
    for l in lines[i + 1:]:
        s = l.strip()
        if not s:
            continue
        parts = s.split()
        if len(parts) < 4 or not parts[0].isalpha() or len(parts[0]) > 2:
            break  # past the block (e.g. 'End ...' or 'JOB DONE')
        # skip if_pos flags (parts[4:7] if present)
        syms.append(parts[0])
        frac.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return syms, np.array(frac)


def load_seed(json_path, system):
    with open(json_path) as f:
        d = json.load(f)
    info = d[system]
    return np.array(info["I_frac"]), info["nI"], info["source"]


def build(system, base="."):
    relax_in = os.path.join(base, system, "activate_empty", "relax.in")
    relax_out = os.path.join(base, system, "activate_empty", "relax.out")
    if not os.path.exists(relax_out):
        raise FileNotFoundError(
            f"{relax_out} not found -- run stage-1 activation first "
            f"(build_activate.py + pw.x).")
    # cell is fixed during activation -> read from the input (pw.x does not echo
    # a CELL_PARAMETERS card for fixed-cell runs)
    cell = read_cell_from_input(relax_in)
    fw_syms, fw_frac = read_activated_framework(relax_out)
    metal = detect_metal(fw_syms)

    seed_path = os.path.join(base, "i2_seed_coords.json")
    i_frac, nI, src = load_seed(seed_path, system)

    outdir = os.path.join(base, system, "activated_i2")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "loaded.in")

    nat = len(fw_syms) + nI
    uniq = species_table(list(fw_syms) + ["I"])
    ntyp = len(uniq)
    degauss = DEGAUSS[system]

    L = []
    L += [
        "&CONTROL",
        "  calculation = 'relax'",
        f"  prefix = '{system}_activated_i2_loaded'",
        "  outdir = './tmp'",
        f"  pseudo_dir = '{PSEUDO_DIR}'",
        "  tprnfor = .true.",
        "  tstress = .true.",
        "  etot_conv_thr = 1.0d-5",
        "  forc_conv_thr = 1.0d-3",
        "/",
        "&SYSTEM",
        "  ibrav = 0",
        f"  nat = {nat}",
        f"  ntyp = {ntyp}",
        f"  ecutwfc = {ECUTWFC}",
        f"  ecutrho = {ECUTRHO}",
        "  occupations = 'smearing'",
        "  smearing = 'mv'",
        f"  degauss = {degauss}",
        "  nspin = 1",
        "  vdw_corr = 'grimme-d3'",
        "  dftd3_version = 4",
        "/",
        "&ELECTRONS",
        "  conv_thr = 1.0d-5",
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
    # framework atoms FROZEN (if_pos 0 0 0), in relax.out order
    for el, f in zip(fw_syms, fw_frac):
        L.append(f"{el:<3s} {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}  0 0 0")
    # I atoms MOVABLE (if_pos 1 1 1), appended last
    for f in i_frac:
        L.append(f"I    {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}  1 1 1")
    L += ["K_POINTS automatic", "  1 1 1 0 0 0", ""]

    with open(path, "w") as fh:
        fh.write("\n".join(L))
    print(f"[{system}] degauss={degauss} nat={nat} "
          f"(fw {len(fw_syms)} frozen + {nI} I movable, seed={src}) -> {path}")
    return path


def main():
    base = "/workspace/dimitry/calc/qe/i2_mof"
    args = sys.argv[1:]
    if not args:
        print(__doc__); sys.exit(1)
    if args[0] == "--check":
        for s in SYSTEMS:
            p = os.path.join(base, s, "activate_empty", "relax.out")
            ok = os.path.exists(p)
            # also check it has a final geometry
            done = False
            if ok:
                with open(p) as fh:
                    done = "JOB DONE" in fh.read()
            print(f"  {s}: relax.out={'YES' if ok else 'NO '}  JOB_DONE={'YES' if done else 'no '}")
        return
    if args[0] == "--all":
        missing = []
        for s in SYSTEMS:
            try:
                build(s, base)
            except FileNotFoundError as e:
                missing.append(s)
                print(f"[{s}] SKIP: {e}")
        if missing:
            print(f"\nNOTE: {len(missing)} systems need stage-1 activation first: {missing}")
        return
    for s in args:
        build(s, base)


if __name__ == "__main__":
    main()
