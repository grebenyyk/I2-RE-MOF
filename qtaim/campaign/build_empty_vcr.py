#!/usr/bin/env python3
"""Build empty-framework RELAX inputs at vc-relax cells + empty_fixed_03 SCF inputs.
For proper binding correction: empty must find its OWN equilibrium at the loaded vc-relax cell."""
import re, os

BASE = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
PP = {"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
      "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF",
      "Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF"}
MASS = {"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"Tb":158.925,"Gd":157.25}

for sys in ["Gd_act","Y_actc"]:
    d = f"{BASE}/{sys}"
    # read vc-relax final geometry
    lines = open(f"{d}/vc-relax.out", errors="ignore").read().splitlines()
    ci = [i for i,l in enumerate(lines) if "CELL_PARAMETERS" in l][-1]
    pi = [i for i,l in enumerate(lines) if "ATOMIC_POSITIONS" in l][-1]
    cell_lines = [lines[ci+1+j].strip() for j in range(3)]
    atoms = []
    for l in lines[pi+1:]:
        m = re.match(r"\s*([A-Z][a-z]?)\s+(-?\d\.\d+)\s+(-?\d\.\d+)\s+(-?\d\.\d+)", l)
        if m and len(atoms) < 80:
            atoms.append((m.group(1), [float(m.group(2)), float(m.group(3)), float(m.group(4))]))
        elif len(atoms) >= 80:
            break
    fw = [(s,f) for s,f in atoms if s != "I"]
    metal = [s for s,_ in atoms if s in ("Y","Tb","Gd")][0]
    uniq = []
    for s,_ in fw:
        if s not in uniq: uniq.append(s)

    # 1. empty RELAX at vc-relax cell
    L = ["&CONTROL",
         "  calculation = 'relax'",
         f"  prefix = '{sys}_empty_vcr'",
         "  outdir = './tmp_evcr'",
         "  pseudo_dir = '/workspace/qe-7.5/pseudo'",
         "  tstress = .true.",
         "  tprnfor = .true.",
         "/",
         "&SYSTEM",
         "  ibrav = 0",
         f"  nat = {len(fw)}",
         f"  ntyp = {len(uniq)}",
         "  ecutwfc = 60",
         "  ecutrho = 800",
         "  occupations = 'smearing'",
         "  smearing = 'mv'",
         "  degauss = 0.03",
         "  nspin = 1",
         "  vdw_corr = 'grimme-d3'",
         "  dftd3_version = 4",
         "/",
         "&ELECTRONS",
         "  conv_thr = 1.0d-8",
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
         "ATOMIC_SPECIES"]
    for s in uniq:
        L.append(f"{s:<3} {MASS[s]:>10.3f}  {PP[s]}")
    L.append("CELL_PARAMETERS angstrom")
    for c in cell_lines:
        L.append(f"  {c}")
    L.append("ATOMIC_POSITIONS crystal")
    for s,f in fw:
        L.append(f"{s:<3} {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}")
    L += ["K_POINTS gamma", ""]
    open(f"{d}/scf_empty_vcr.in", "w").write("\n".join(L))
    print(f"{sys}: empty RELAX at vc-relax cell ({len(fw)} fw atoms, metal={metal})")

    # 2. empty_fixed_03 SCF
    s2 = open(f"{d}/scf_empty.in").read()
    s2 = s2.replace("degauss = 0.01", "degauss = 0.03")
    s2 = s2.replace("prefix = '", "prefix = 'fix03_")
    open(f"{d}/scf_empty_fix03.in", "w").write(s2)
    print(f"{sys}: empty_fixed_03 SCF built (degauss 0.03)")
