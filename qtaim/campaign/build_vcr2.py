#!/usr/bin/env python3
"""Build empty_vcr2 RELAX inputs (damped dynamics) from step-49 BFGS geometry.
Gd/Y: extract last ATOMIC_POSITIONS from vcr.out.
Tb: cell-swap of fix03 (same approach that worked for SCF).
All use ion_dynamics='damp' to avoid BFGS stalling."""
import re

BASE = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
PP = {"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
      "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF",
      "Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF"}
MASS = {"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"Tb":158.925,"Gd":157.25}

def extract_last_geometry(outpath):
    """Extract last ATOMIC_POSITIONS block from QE relax output."""
    lines = open(outpath, errors="ignore").read().splitlines()
    pos_indices = [i for i,l in enumerate(lines) if "ATOMIC_POSITIONS" in l]
    if not pos_indices:
        return None
    pi = pos_indices[-1]
    atoms = []
    for l in lines[pi+1:]:
        m = re.match(r"\s*([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d\.\d+)", l)
        if m and len(atoms) < 72:
            atoms.append((m.group(1), [float(m.group(2)), float(m.group(3)), float(m.group(4))]))
        elif len(atoms) >= 72:
            break
    return atoms

def get_cell_from_input(inpath):
    """Extract CELL_PARAMETERS from an existing input file."""
    s = open(inpath).read()
    cell_block = re.search(r"CELL_PARAMETERS\s+angstrom\n((?:\s*-?[\d.]+.*\n){3})", s)
    if cell_block:
        return [line.strip() for line in cell_block.group(1).strip().split("\n")]
    return None

def write_relax_input(sys, atoms, cell_lines, metal):
    uniq = []
    for s,_ in atoms:
        if s not in uniq: uniq.append(s)
    L = [
        "&CONTROL",
        "  calculation = 'relax'",
        f"  prefix = '{sys}_empty_vcr2'",
        "  outdir = './tmp_evcr2'",
        "  pseudo_dir = '/workspace/qe-7.5/pseudo'",
        "  tstress = .true.",
        "  tprnfor = .true.",
        "/",
        "&SYSTEM",
        "  ibrav = 0",
        f"  nat = {len(atoms)}",
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
        "  nstep = 200",
        "/",
        "&IONS",
        "  ion_dynamics = 'damp'",
        "  forc_conv_thr = 1.0d-3",
        "/",
        "ATOMIC_SPECIES",
    ]
    for s in uniq:
        L.append(f"{s:<3} {MASS[s]:>10.3f}  {PP[s]}")
    L.append("CELL_PARAMETERS angstrom")
    for c in cell_lines:
        L.append(f"  {c}")
    L.append("ATOMIC_POSITIONS crystal")
    for s,f in atoms:
        L.append(f"{s:<3} {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}")
    L += ["K_POINTS gamma", ""]
    outpath = f"{BASE}/{sys}/scf_empty_vcr2.in"
    open(outpath, "w").write("\n".join(L))
    return len(atoms)

for sys in ["Gd_act", "Y_actc", "Tb_act"]:
    metal = {"Gd_act":"Gd","Y_actc":"Y","Tb_act":"Tb"}[sys]

    if sys in ("Gd_act", "Y_actc"):
        # Extract last geometry from BFGS output
        atoms = extract_last_geometry(f"{BASE}/{sys}/scf_empty_vcr.out")
        cell_lines = get_cell_from_input(f"{BASE}/{sys}/scf_empty_vcr.in")
        if atoms is None or cell_lines is None:
            print(f"{sys}: FAILED to extract geometry")
            continue
    else:
        # Tb: use fix03 input atoms + vc-relax cell
        fix03 = open(f"{BASE}/{sys}/scf_empty_fix03.in").read()
        cell_lines = get_cell_from_input(f"{BASE}/{sys}/scf_empty_vcr.in")
        if cell_lines is None:
            print(f"{sys}: FAILED to extract cell")
            continue
        # parse atoms from fix03
        lines_fix = fix03.splitlines()
        pi = [i for i,l in enumerate(lines_fix) if "ATOMIC_POSITIONS" in l][-1]
        atoms = []
        for l in lines_fix[pi+1:]:
            m = re.match(r"\s*([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d\.\d+)\s+(-?\d\.\d+)", l)
            if m and len(atoms) < 72:
                atoms.append((m.group(1), [float(m.group(2)), float(m.group(3)), float(m.group(4))]))
            elif len(atoms) >= 72:
                break

    nat = write_relax_input(sys, atoms, cell_lines, metal)
    print(f"{sys}: vcr2 RELAX input built ({nat} fw atoms, damp dynamics)")
