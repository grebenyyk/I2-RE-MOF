#!/usr/bin/env python3
"""Generate full-vc-relax inputs for the EMPTY frameworks (no I2) of the 3 study
MOFs x 3 starting geometries (st/act/vc = str/act/deep). Goal: let both atoms
AND cell relax freely (cell_dofree='all') from each empty starting geometry, to
see where the guest-free framework wants to go.

Source geometry: each variant's scf_empty.in (framework only, nat=72). We rebuild
the namelists (vc-relax + manuscript settings) and keep the ATOMIC_SPECIES /
CELL_PARAMETERS / ATOMIC_POSITIONS cards verbatim. Writes <variant>/vc_empty.in
(never overwrites; never touches scf_empty.in or the loaded vc-relax.in).

Settings (manuscript geom-opt): ecutwfc=50, ecutrho=400, cold smearing (mv)
degauss=0.01, conv_thr=1e-6; vc-relax etot_conv=1e-5, forc_conv=1e-3,
press_conv=0.5; Gamma. nspin=1, f-in-PP (no Hubbard), PBE+D3(BJ) v4.

NOTE: the primary campaign geometry relaxations (*/loaded.in, */activate_empty/
relax.in, full_relax/*, all per-metal relaxes) used ecutwfc=50, ecutrho=400 --
matching the manuscript geom-opt settings. The 60/800 files are a DIFFERENT
category (detailed SCF scf_empty.in for QTAIM charge density, the high-precision
loaded-cell refine qtaim/campaign/vc-relax.in, scf_empty_vcr2.in), NOT the
geom-opt setting. So 50/400 here is fully consistent with the campaign geom-opts.
The starting GEOMETRY is read from scf_empty.in (the converged empty framework)
and is fine regardless of the ecut it was originally SCF'd at.
"""
import os, re, numpy as np

BASE = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
VARIANTS = ["Y_str","Y_actc","Y_deep","Gd_str","Gd_act","Gd_deep",
            "Tb_str","Tb_act","Tb_deep"]

# ---- settings (manuscript geometry-optimization) ----
ECUTWFC, ECUTRHO = 50, 400        # Ry ; campaign files used 60/800 (see NOTE)
DEGAUSS = 0.01                     # cold smearing
CONV    = "1.0d-6"                 # SCF conv_thr
ETOT    = "1.0d-5"                 # vc-relax energy threshold
FORC    = "1.0d-3"                 # vc-relax force threshold
PRESS   = 0.5                      # kbar, cell pressure threshold
MIXBETA = 0.2                      # empty framework is closed-shell; scf_empty validated at 0.2

def cards_from_scf_empty(f):
    """Return verbatim ATOMIC_SPECIES + CELL_PARAMETERS + ATOMIC_POSITIONS cards."""
    txt = open(f).read()
    i = txt.find("ATOMIC_SPECIES")
    assert i != -1, f"no ATOMIC_SPECIES in {f}"
    cards = txt[i:].split("K_POINTS")[0].rstrip()   # drop the scf K_POINTS card
    return cards

def cell_volume(cards):
    m = re.search(r"CELL_PARAMETERS.*?\n((?:\s*-?\d[\d.\s-]*\n){3})", cards, re.DOTALL)
    if not m: return float("nan")
    c = np.array([[float(x) for x in ln.split()] for ln in m.group(1).splitlines()])
    return abs(np.linalg.det(c))

for v in VARIANTS:
    src = f"{BASE}/{v}/scf_empty.in"
    out = f"{BASE}/{v}/vc_empty.in"
    if not os.path.exists(src):
        print(f"  SKIP {v}: no scf_empty.in"); continue
    if os.path.exists(out):
        print(f"  SKIP {v}: vc_empty.in already exists (not overwriting)"); continue
    cards = cards_from_scf_empty(src)
    vol = cell_volume(cards)
    nat  = len(re.findall(r"^\s*[A-Z][a-z]?\s+-?\d", cards, re.MULTILINE))
    prefix = f"{v}_empty_vc"
    nl = (
f"""&CONTROL
  calculation = 'vc-relax'
  prefix = '{prefix}'
  outdir = './tmp_empty_vc'
  pseudo_dir = '/workspace/qe-7.5/pseudo'
  tprnfor = .true.
  tstress = .true.
  etot_conv_thr = {ETOT}
  forc_conv_thr = {FORC}
/
&SYSTEM
  ibrav = 0
  nat = 72
  ntyp = 4
  ecutwfc = {ECUTWFC}
  ecutrho = {ECUTRHO}
  occupations = 'smearing'
  smearing = 'mv'
  degauss = {DEGAUSS}
  nspin = 1
  vdw_corr = 'grimme-d3'
  dftd3_version = 4
/
&ELECTRONS
  conv_thr = {CONV}
  diagonalization = 'david'
  diago_david_ndim = 2
  mixing_mode = 'local-TF'
  mixing_beta = {MIXBETA}
  mixing_ndim = 16
  electron_maxstep = 300
/
&IONS
  ion_dynamics = 'bfgs'
/
&CELL
  cell_dynamics = 'bfgs'
  cell_dofree = 'all'
  press_conv_thr = {PRESS}
  cell_factor = 2.0
/
""")
    open(out, "w").write(nl + cards + "\nK_POINTS gamma\n")
    print(f"  wrote {v}/vc_empty.in   (nat={nat}, V0={vol:.1f} A^3, ecut={ECUTWFC}/{ECUTRHO})")
print("done.")
