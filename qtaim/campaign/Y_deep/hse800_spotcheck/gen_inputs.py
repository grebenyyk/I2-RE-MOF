#!/usr/bin/env python3
"""Generate HSE06 inputs at the production PBE-density settings (60/800/Gamma/1e-8)
on the IDENTICAL vc geometry as Y_deep/rho17.cube (cell 8.7684... from hse/scf.in).
Emits loaded (nat=80) and empty-host (nat=72, I removed) for Bader charge-transfer."""
import re
SRC="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/hse/scf.in"
OUT="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/hse800_spotcheck"
txt=open(SRC).read()
# parse CELL + ATOMIC_POSITIONS
cell=txt.split("CELL_PARAMETERS")[1].splitlines()[1:4]
apos=[]
for l in txt.split("ATOMIC_POSITIONS")[1].split("K_POINTS")[0].splitlines():
    m=re.match(r"\s*([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",l)
    if m: apos.append((m.group(1),m.group(2),m.group(3),m.group(4)))
assert len(apos)==80, len(apos)
# ATOMIC_SPECIES block (drop I for empty)
SPEC_L="C       12.011  C.pbe-n-kjpaw_psl.1.0.0.UPF\nH        1.008  H.pbe-kjpaw.UPF\nO       15.999  O.pbe-kjpaw.UPF\nY       88.906  Y.pbe-spn-kjpaw_psl.1.0.0.UPF"
SPEC_F=SPEC_L+"\nI      126.904  I.pbe-n-kjpaw_psl.1.0.0.UPF"
def namelist(nat,ntyp,spec,prefix,outdir):
    return f"""&CONTROL
  calculation = 'scf'
  prefix = '{prefix}'
  outdir = '{outdir}'
  pseudo_dir = '/workspace/qe-7.5/pseudo'
  verbosity = 'high'
  wf_collect = .true.
  disk_io = 'high'
/
&SYSTEM
  ibrav = 0
  nat = {nat}
  ntyp = {ntyp}
  ecutwfc = 60
  ecutrho = 800
  occupations = 'smearing'
  smearing = 'mv'
  degauss = 0.01
  nspin = 1
  input_dft = 'hse'
  exxdiv_treatment = 'gygi-baldereschi'
  x_gamma_extrapolation = .true.
  exx_fraction = 0.25
  screening_parameter = 0.106
  ecutfock = 100
/
&ELECTRONS
  conv_thr = 1.0d-8
  diagonalization = 'david'
  diago_david_ndim = 2
  mixing_mode = 'local-TF'
  mixing_beta = 0.05
  mixing_ndim = 24
  electron_maxstep = 300
/
ATOMIC_SPECIES
{spec}
CELL_PARAMETERS angstrom
{chr(10).join(cell)}
ATOMIC_POSITIONS crystal
"""
def write_atoms(f,atoms):
    for s,x,y,z in atoms: f.write(f"  {s:<2s} {x} {y} {z}\n")
# loaded: all 80
with open(f"{OUT}/hse800_vc_loaded.in","w") as f:
    f.write(namelist(80,5,SPEC_F,"hse800_vc","./tmp_loaded"))
    write_atoms(f,apos)
    f.write("K_POINTS automatic\n  1 1 1 0 0 0\n")
# empty: drop 8 I -> 72 atoms (C/H/O/Y), ntyp=4
empty=[a for a in apos if a[0]!="I"]
assert len(empty)==72, len(empty)
with open(f"{OUT}/hse800_vc_empty.in","w") as f:
    f.write(namelist(72,4,SPEC_L,"hse800_vc_empty","./tmp_empty"))
    write_atoms(f,empty)
    f.write("K_POINTS automatic\n  1 1 1 0 0 0\n")
print(f"wrote loaded (nat=80) and empty (nat=72)")
print(f"cell[0]={cell[0].strip()}")
print(f"empty species: {sorted(set(a[0] for a in empty))}")
