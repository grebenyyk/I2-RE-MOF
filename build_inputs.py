#!/usr/bin/env python3
"""
Build Quantum ESPRESSO inputs for I2@MOF interaction-energy calculations.

For one I2-loaded MOF unit cell (a POSCAR from poscars_unitcell/), produces:
  <outdir>/empty_mof.in  : SCF on the bare framework (I atoms removed), framework frozen
  <outdir>/loaded.in     : relax of the full I2@MOF, framework frozen, only I atoms move

Interaction energy (per MOF):
  E_int   = E(loaded) - E(empty) - N_I2 * E(I2_vac)
  E_per_I2 = E_int / N_I2
with the shared I2 vacuum reference E(I2_vac) = -757.80004425 Ry (PBE+D3(BJ), 50/400,
I.pbe-n-kjpaw_psl; computed once in i2_ce-uio-66/01_i2_vac).

Setup mirrors the converged i2_ce-uio-66 loaded template (smearing mv degauss 0.01,
diagonalization paro, local-TF mixing beta 0.2 ndim 16, conv_thr 1d-5, 50/400, D3(BJ)
dftd3_version=4, Gamma point, framework frozen via if_pos 0 0 0, I atoms if_pos 1 1 1).

Metal treatment (the only physics choice that differs from the Ce4+ reference, which was
4f0 / nspin=1 / no-U):
  - closed shell  (Cd2+ 4d10, Y3+ 4d0):               nspin=1, no Hubbard
  - open shell f  (Tb3+ 4f8, Gd3+ 4f7, U4+ 5f2):       nspin=2, Hubbard U on f, starting_magnetization
The f-electron energy largely cancels in E_int PROVIDED the metal f-occupation is held
identical between empty and loaded (same magnetic state, same U) -- I2 adsorption does not
redox the metal. Empty and loaded are therefore generated with identical nspin/U/mag.

Usage:
  python build_inputs.py <POSCAR> <outdir> <label>
  python build_inputs.py --pilot        # builds the Y (13) and Gd (08) pilot systems
"""
import os, sys
import numpy as np
import ase.io

PSEUDO_DIR = "/workspace/qe-7.5/pseudo"
ECUTWFC = int(os.environ.get("ECUTWFC", 50))
ECUTRHO = int(os.environ.get("ECUTRHO", 400))
# Smearing/mixing are env-overridable. Default 0.01/0.2 works for most systems, but some
# I2@MOF loaded relaxes hit a geometry where the gap < degauss -> Fermi minimization fails
# -> NaN -> GPU bus error (SIGABRT). For those, rerun with DEGAUSS=0.005 MIXING_BETA=0.15.
DEGAUSS = os.environ.get("DEGAUSS", "0.01")
MIXING_BETA = os.environ.get("MIXING_BETA", "0.2")

# --- pseudopotentials + masses (PBE PAW, PSlibrary 1.0.0 legacy table) ---
COMMON = {
    "C": ("C.pbe-n-kjpaw_psl.1.0.0.UPF", 12.011),
    "H": ("H.pbe-kjpaw.UPF", 1.008),
    "O": ("O.pbe-kjpaw.UPF", 15.999),
    "N": ("N.pbe-n-kjpaw_psl.1.0.0.UPF", 14.007),
    "I": ("I.pbe-n-kjpaw_psl.1.0.0.UPF", 126.904),
}
# Metal-specific: (pseudo, mass, f_in_core, nspin, hubbard_orbital, hubbard_U_eV, starting_magnetization)
#
# IMPORTANT: the PSlibrary "spdn" PPs for the lanthanides/actinides (Tb, Gd, U) are
# f-IN-CORE: the 4f/5f shell is frozen into a non-spin-polarized core and is NOT a
# variational valence channel (QE lists only 5S6S5P6P5D as manifolds -> no f to apply U
# on, no f moment). With such a PP the trivalent ion's VALENCE is closed shell
# (Ln3+: 5d0 6s0 6p0; the 3 valence e- transfer to O/linkers), so the correct treatment
# is nspin=1, NO Hubbard U -- exactly like Ce4+ 4f0 in i2_ce-uio-66. The frozen f-shell is
# chemically inert and identical in empty vs loaded, so it cancels exactly in E_int.
#
# f_in_core=True  -> closed-shell valence: nspin=1, no U, no starting_magnetization
# f_in_core=False -> f-in-VALENCE PP (NOT what PSlibrary ships): nspin=2 + U on f + mag.
#                    The U/mag values below are literature starting points kept for that case.
METALS = {
    "Y":  ("Y.pbe-spn-kjpaw_psl.1.0.0.UPF",  88.906, True,  1, None,  None, None),   # Y3+ 4d0 closed
    "Cd": ("Cd.pbe-n-kjpaw_psl.1.0.0.UPF",  112.414, True,  1, None,  None, None),   # Cd2+ 4d10 closed (no spn variant exists)
    "Gd": ("Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF", 157.25, True,  1, None,  None, None),   # Gd3+ 4f7 frozen-in-core
    "Tb": ("Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF", 158.925, True,  1, None,  None, None),  # Tb3+ 4f8 frozen-in-core
    "U":  ("U.pbe-spdn-kjpaw_psl.1.0.0.UPF",  238.029, True,  1, None,  None, None),   # U4+ 5f2 frozen-in-core (VERIFY f-count vs oxidation)
}
# If you later switch to f-in-VALENCE PPs, use these:
METALS_FVALENCE = {
    "Gd": ("Gd-4f-in-valence.UPF", 157.25, False, 2, "4f", 6.0, 0.6),   # Gd3+ 4f7, ~7 muB
    "Tb": ("Tb-4f-in-valence.UPF", 158.925, False, 2, "4f", 6.0, 0.5),  # Tb3+ 4f8, ~6 muB
    "U":  ("U-5f-in-valence.UPF",  238.029, False, 2, "5f", 4.0, 0.2),  # U4+ 5f2, ~2 muB
}
FRAMEWORK_ELS = set(COMMON) - {"I"}   # C, H, O, N


def detect_metal(symbols):
    metals = [s for s in symbols if s not in COMMON]
    assert len(metals) > 0, "no metal found (only C/H/O/N/I present)"
    assert len(set(metals)) == 1, f"multiple metals {set(metals)} -- mono-metallic only"
    return metals[0]


def species_table(symbols):
    """Return ordered unique species list (metal first for stable HUBBARD indexing)."""
    uniq = []
    for s in symbols:
        if s not in uniq:
            uniq.append(s)
    return uniq


def pseudo_mass(el, metal):
    if el == metal:
        m = METALS[metal]
        return m[0], m[1]
    return COMMON[el]


def build(poscar, outdir, label):
    atoms = ase.io.read(poscar)
    cell = np.array(atoms.cell)
    frac = atoms.get_scaled_positions(wrap=False)
    syms = list(atoms.get_chemical_symbols())
    metal = detect_metal(syms)
    f_in_core = METALS[metal][2]
    nspin, orb, Uval, startmag = METALS[metal][3:7]
    n_I = syms.count("I")
    assert n_I % 2 == 0, f"odd number of I atoms ({n_I})"
    N_I2 = n_I // 2

    os.makedirs(outdir, exist_ok=True)

    # ---- empty framework: drop I atoms, keep order ----
    fw_idx = [i for i, s in enumerate(syms) if s != "I"]
    fw_syms = [syms[i] for i in fw_idx]
    fw_frac = frac[fw_idx]

    # ---- loaded: all atoms, I movable ----
    ld_syms = syms
    ld_frac = frac

    def write_input(path, calc, syms_list, frac_list, movable_test, prefix):
        uniq = species_table(syms_list)
        nat = len(syms_list)
        ntyp = len(uniq)
        metal_typ_idx = uniq.index(metal) + 1   # 1-based for starting_magnetization / HUBBARD
        lines = []
        lines += [
            "&CONTROL",
            f"  calculation = '{calc}'",
            f"  prefix = '{prefix}'",
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
            f"  degauss = {DEGAUSS}",
            f"  nspin = {nspin}",
            "  vdw_corr = 'grimme-d3'",
            "  dftd3_version = 4",
        ]
        if nspin == 2 and startmag is not None:
            lines.append(f"  starting_magnetization({metal_typ_idx}) = {startmag}")
        lines += [
            "/",
            "&ELECTRONS",
            "  conv_thr = 1.0d-5",
            "  diagonalization = 'david'",
            "  diago_david_ndim = 4",
            "  mixing_mode = 'local-TF'",
            f"  mixing_beta = {MIXING_BETA}",
            "  mixing_ndim = 16",
            "  electron_maxstep = 300",
            "/",
        ]
        if calc == "relax":
            lines += ["&IONS", "  ion_dynamics = 'bfgs'", "/"]
        # ATOMIC_SPECIES
        lines += ["ATOMIC_SPECIES"]
        for el in uniq:
            p, m = pseudo_mass(el, metal)
            lines.append(f"{el:<3s} {m:>10.3f}  {p}")
        # CELL
        lines += ["CELL_PARAMETERS angstrom"]
        for r in cell:
            lines.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}")
        # POSITIONS (crystal/fractional); add if_pos for relax
        lines += ["ATOMIC_POSITIONS crystal"]
        for el, f in zip(syms_list, frac_list):
            row = f"{el:<3s} {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}"
            if calc == "relax":
                row += "  " + ("1 1 1" if movable_test(el) else "0 0 0")
            lines.append(row)
        # HUBBARD card
        if orb is not None and Uval is not None:
            lines += [f"HUBBARD {{atomic}}", f" U {metal}-{orb} {Uval}"]
        lines += ["K_POINTS automatic", "  1 1 1 0 0 0", ""]
        with open(path, "w") as fh:
            fh.write("\n".join(lines))
        return path

    ep = write_input(os.path.join(outdir, "empty_mof.in"), "scf",
                     fw_syms, fw_frac, lambda el: False, f"{label}_empty")
    lp = write_input(os.path.join(outdir, "loaded.in"), "relax",
                     ld_syms, ld_frac, lambda el: el == "I", f"{label}_loaded")
    print(f"[{label}] metal={metal} nspin={nspin} U={Uval} N_I2={N_I2}  ->  {ep}, {lp}")
    return dict(metal=metal, nspin=nspin, U=Uval, N_I2=N_I2,
                empty=ep, loaded=lp, outdir=outdir, label=label)


if __name__ == "__main__":
    base = "/workspace/dimitry/calc/qe/i2_mof"
    pos = os.path.join(base, "poscars_unitcell")
    # label -> POSCAR filename. U systems (16, 18) skipped per user (oxidation/f-count risk).
    MOFS = {
        "13_Y":  "POSCAR_13_2020_Y_pts_3_ASR_1_I2loaded_1000Pa_UC",
        "08_Gd": "POSCAR_08_2020_Gd_pts_3_ASR_1_I2loaded_1000Pa_UC",
        "05_Cd": "POSCAR_05_oc8b00378_si_004_ASR_pacman_I2loaded_1000Pa_UC",
        "06_Tb": "POSCAR_06_2019_Tb_pts_3_ASR_1_I2loaded_1000Pa_UC",
        "07_Cd": "POSCAR_07_oc8b00378_si_006_ASR_pacman_I2loaded_1000Pa_UC",
        "10_Cd": "POSCAR_10_oc8b00378_si_003_ASR_pacman_I2loaded_1000Pa_UC",
        "12_Tb": "POSCAR_12_2020_Tb_nan_3_ASR_7_I2loaded_1000Pa_UC",
        "17_Cd": "POSCAR_17_2016_Cd_nan_3_ASR_1_I2loaded_1000Pa_UC",
    }
    if "--pilot" in sys.argv:
        for sub in ("13_Y", "08_Gd"):
            build(os.path.join(pos, MOFS[sub]), os.path.join(base, sub), sub)
    elif "--campaign" in sys.argv:
        # build the 6 remaining MOFs (pilot 13_Y/08_Gd already done; U 16/18 skipped)
        for sub in ("05_Cd", "06_Tb", "07_Cd", "10_Cd", "12_Tb", "17_Cd"):
            build(os.path.join(pos, MOFS[sub]), os.path.join(base, sub), sub)
    else:
        if len(sys.argv) != 4:
            print(__doc__); sys.exit(1)
        build(sys.argv[1], sys.argv[2], sys.argv[3])
