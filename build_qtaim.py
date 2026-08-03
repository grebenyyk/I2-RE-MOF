#!/usr/bin/env python3
"""Build tight single-point SCF + pp.x (all-electron density, plot_num=21) inputs for QTAIM,
starting from a system's POSCAR_final (relaxed geometry).

Usage: python build_qtaim.py <system> <ecutwfc> <ecutrho> <k1 k2 k3> <outdir>
e.g.  python build_qtaim.py 13_Y 60 600 "2 2 2" qtaim/13_Y_rho600_2x2x2
"""
import sys, os
import numpy as np, ase.io

PSEUDO_DIR = "/workspace/qe-7.5/pseudo"
DEGAUSS = {"13_Y": 0.01, "06_Tb": 0.01, "08_Gd": 0.01}  # match each system's relax
PSEUDOS = {"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
           "N":"N.pbe-n-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF",
           "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF",
           "Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF","Cd":"Cd.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS = {"C":12.011,"H":1.008,"O":15.999,"N":14.007,"I":126.904,"Y":88.906,
        "Gd":157.25,"Tb":158.925,"Cd":112.414}


def build(system, ecutwfc, ecutrho, ks, outdir):
    # prefer the uniform_chain optimal geometry; fall back to system POSCAR_final
    geo = f"{system}/uniform_chain/POSCAR_final"
    if not os.path.exists(geo):
        geo = f"{system}/POSCAR_final"
    atoms = ase.io.read(geo)
    cell = np.array(atoms.cell)
    frac = atoms.get_scaled_positions(wrap=False)
    syms = list(atoms.get_chemical_symbols())
    uniq = []
    for s in syms:
        if s not in uniq:
            uniq.append(s)
    deg = DEGAUSS.get(system, 0.01)
    os.makedirs(outdir, exist_ok=True)
    pref = os.path.basename(outdir.rstrip("/"))
    L = ["&CONTROL", "  calculation = 'scf'", f"  prefix = '{pref}'",
         "  outdir = './tmp'", f"  pseudo_dir = '{PSEUDO_DIR}'",
         "  tprnfor = .true.", "  tstress = .true.", "/",
         "&SYSTEM", "  ibrav = 0", f"  nat = {len(syms)}", f"  ntyp = {len(uniq)}",
         f"  ecutwfc = {ecutwfc}", f"  ecutrho = {ecutrho}",
         "  occupations = 'smearing'", "  smearing = 'mv'", f"  degauss = {deg}",
         "  nspin = 1", "  vdw_corr = 'grimme-d3'", "  dftd3_version = 4", "/",
         "&ELECTRONS", "  conv_thr = 1.0d-8", "  diagonalization = 'david'",
         "  diago_david_ndim = 4", "  mixing_mode = 'local-TF'", "  mixing_beta = 0.2",
         "  mixing_ndim = 16", "  electron_maxstep = 300", "/",
         "ATOMIC_SPECIES"]
    for s in uniq:
        L.append(f"{s:<3s} {MASS[s]:>10.3f}  {PSEUDOS[s]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell:
        L.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for s, f in zip(syms, frac):
        L.append(f"{s:<3s} {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}")
    L += ["K_POINTS automatic", f"  {ks[0]} {ks[1]} {ks[2]} 0 0 0", ""]
    scf_in = os.path.join(outdir, "scf.in")
    open(scf_in, "w").write("\n".join(L))

    # pp.x input: all-electron charge (plot_num=21) -> Gaussian cube (output_format=6).
    # Two namelists: &INPUTPP selects the quantity, &PLOT writes the cube.
    pp = ["&INPUTPP", f"  prefix = '{pref}'", "  outdir = './tmp'",
          "  plot_num = 21", "  filplot = 'ae_rho'", "/",
          "&PLOT", "  nfile = 1", "  filepp(1) = 'ae_rho'", "  weight(1) = 1.0",
          "  iflag = 3", "  output_format = 6", "  fileout = 'ae_rho.cube'", "/"]
    open(os.path.join(outdir, "pp.in"), "w").write("\n".join(pp))
    print(f"{pref}: scf.in (ecutwfc={ecutwfc} ecutrho={ecutrho} k={ks}) + pp.in (plot_num=21 cube)")


if __name__ == "__main__":
    s, ew, er, ks_str, out = sys.argv[1:6]
    build(s, int(ew), int(er), [int(x) for x in ks_str.split()], out)
