#!/usr/bin/env python3
"""Build continued-BFGS (tight) relax inputs for the 5 strongest binders, starting from
the relaxed POSCAR_final of each system. Settings:
  - nstep=200            (previous runs hit the default 50-step cap -> the real bottleneck)
  - conv_thr=1e-6        (tighter SCF -> clean forces for BFGS force criterion)
  - etot_conv_thr=1e-5, forc_conv_thr=1e-3  (standard "converged" criteria)
  - david + diago_david_ndim=4, local-TF mixing
  - degauss per system MATCHING its empty reference (so ΔE stays valid):
      13_Y/06_Tb/08_Gd -> 0.01 (mixing 0.2),  12_Tb/17_Cd -> 0.005 (mixing 0.15)
  - framework frozen (if_pos 0), only I atoms move (if_pos 1)
  - start coords = POSCAR_final (the last relaxed geometry)
"""
import numpy as np, ase.io

PSEUDO_DIR = "/workspace/qe-7.5/pseudo"
SETTINGS = {  # label: (degauss, mixing_beta)
    "13_Y":  (0.01,  0.20),
    "06_Tb": (0.01,  0.20),
    "08_Gd": (0.01,  0.20),
    "12_Tb": (0.005, 0.15),
    "17_Cd": (0.005, 0.15),
}
PSEUDOS = {
    "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF", "H": "H.pbe-kjpaw.UPF", "O": "O.pbe-kjpaw.UPF",
    "N": "N.pbe-n-kjpaw_psl.1.0.0.UPF", "I": "I.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Y": "Y.pbe-spn-kjpaw_psl.1.0.0.UPF", "Gd": "Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "Tb": "Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF", "Cd": "Cd.pbe-n-kjpaw_psl.1.0.0.UPF",
}
MASS = {"C":12.011,"H":1.008,"O":15.999,"N":14.007,"I":126.904,
        "Y":88.906,"Gd":157.25,"Tb":158.925,"Cd":112.414}


def build(m):
    atoms = ase.io.read(f"{m}/POSCAR_final")
    cell = np.array(atoms.cell)
    frac = atoms.get_scaled_positions(wrap=False)
    syms = list(atoms.get_chemical_symbols())
    degauss, mixb = SETTINGS[m]
    uniq = []
    for s in syms:
        if s not in uniq:
            uniq.append(s)
    L = [
        "&CONTROL",
        "  calculation = 'relax'",
        f"  prefix = '{m}_tight'",
        "  outdir = './tmp'",
        f"  pseudo_dir = '{PSEUDO_DIR}'",
        "  tprnfor = .true.", "  tstress = .true.",
        "  nstep = 200",
        "  etot_conv_thr = 1.0d-5", "  forc_conv_thr = 1.0d-3",
        "/",
        "&SYSTEM",
        "  ibrav = 0", f"  nat = {len(syms)}", f"  ntyp = {len(uniq)}",
        "  ecutwfc = 50", "  ecutrho = 400",
        "  occupations = 'smearing'", "  smearing = 'mv'", f"  degauss = {degauss}",
        "  nspin = 1", "  vdw_corr = 'grimme-d3'", "  dftd3_version = 4",
        "/",
        "&ELECTRONS",
        "  conv_thr = 1.0d-6",
        "  diagonalization = 'david'", "  diago_david_ndim = 4",
        "  mixing_mode = 'local-TF'", f"  mixing_beta = {mixb}", "  mixing_ndim = 16",
        "  electron_maxstep = 300",
        "/",
        "&IONS", "  ion_dynamics = 'bfgs'", "/",
        "ATOMIC_SPECIES",
    ]
    for s in uniq:
        L.append(f"{s:<3s} {MASS[s]:>10.3f}  {PSEUDOS[s]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell:
        L.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for s, f in zip(syms, frac):
        L.append(f"{s:<3s} {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}  " + ("1 1 1" if s == "I" else "0 0 0"))
    L += ["K_POINTS automatic", "  1 1 1 0 0 0", ""]
    open(f"{m}/loaded_tight.in", "w").write("\n".join(L))
    print(f"{m}: loaded_tight.in  nat={len(syms)} degauss={degauss} nstep=200 conv_thr=1e-6")


if __name__ == "__main__":
    for m in ["13_Y", "06_Tb", "12_Tb", "08_Gd", "17_Cd"]:
        build(m)
