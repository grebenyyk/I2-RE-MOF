#!/usr/bin/env python3
"""Stage-1 inputs for the 442 Raman campaign: tight I-only re-relaxation of the
four production winner geometries (frozen-cell campaign). Framework atoms get
if_pos 0 0 0; the 8 iodines relax at production smearing 0.02 with conv_thr
1e-9 / forc 1e-5 Ry/bohr (the winners were relaxed at 0.04 / 5e-3 - far too
loose for finite-difference phonons). Geometry seed = the verified winner CIFs.
"""
import math
import numpy as np
import os

CIFS = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results_442/landscape_3variants/cifs'
BASE = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442'
SYS = [('st', 'st_poly_lS5'), ('vc', 'vc_poly_lS3'), ('act', 'act_poly_lS2'), ('actmol', 'act_mol_rd')]
SPECIES = [('C', 12.011, 'C.pbe-n-kjpaw_psl.1.0.0.UPF'),
           ('H', 1.008, 'H.pbe-kjpaw.UPF'),
           ('O', 15.999, 'O.pbe-kjpaw.UPF'),
           ('Y', 88.906, 'Y.pbe-spn-kjpaw_psl.1.0.0.UPF'),
           ('I', 126.904, 'I.pbe-n-kjpaw_psl.1.0.0.UPF')]

TPL = """&CONTROL
  calculation = 'relax'
  prefix = 'raman442_{tag}_rerelax'
  outdir = './tmp_rerelax'
  pseudo_dir = '/workspace/qe-7.5/pseudo'
  tprnfor = .true.
  tstress = .true.
  etot_conv_thr = 1.0d-6
  forc_conv_thr = 1.0d-5
  nstep = 300
/
&SYSTEM
  ibrav = 0
  nat = 80
  ntyp = 5
  ecutwfc = 50
  ecutrho = 400
  occupations = 'smearing'
  smearing = 'mv'
  degauss = 0.02
  nspin = 1
  input_dft = 'pbe'
  vdw_corr = 'grimme-d3'
  dftd3_version = 4
/
&ELECTRONS
  conv_thr = 1.0d-9
  diagonalization = 'david'
  diago_david_ndim = 2
  mixing_mode = 'local-TF'
  mixing_beta = 0.2
  mixing_ndim = 16
  electron_maxstep = 300
/
&IONS
  ion_dynamics = 'bfgs'
/
ATOMIC_SPECIES
""" + ''.join(f'  {e}  {m:.3f}  {p}\n' for e, m, p in SPECIES)

for var, tag in SYS:
    s = open(f'{CIFS}/l442_{tag}_relaxed.cif').read()
    lab = s.index('_atom_site_fract_z')
    par = {}
    for l in s[:lab].splitlines():
        q = l.split()
        if len(q) == 2 and q[0].startswith('_cell_length'):
            par[q[0][-1]] = float(q[1])
        if len(q) == 2 and q[0].startswith('_cell_angle'):
            par[q[0].split('_')[-1]] = float(q[1])
    a, b, c = par['a'], par['b'], par['c']
    al, be, ga = (math.radians(par[k]) for k in ('alpha', 'beta', 'gamma'))
    M = np.array([
        [a, 0, 0],
        [b * math.cos(ga), b * math.sin(ga), 0],
        [c * math.cos(be), c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga),
         c * math.sqrt(1 - math.cos(be)**2 - math.cos(al)**2 - math.cos(ga)**2
                       + 2 * math.cos(be) * math.cos(al) * math.cos(ga)) / math.sin(ga)]])
    fr = []
    for l in s[lab:].splitlines():
        q = l.split()
        if len(q) >= 5 and q[1] in ('C', 'H', 'O', 'Y', 'I'):
            fr.append((q[1], float(q[2]), float(q[3]), float(q[4])))
    assert len(fr) == 80, (tag, len(fr))
    outdir = f'{BASE}/{tag}'
    with open(f'{outdir}/rerelax_{tag}.in', 'w') as f:
        f.write(TPL.format(tag=tag))
        f.write('CELL_PARAMETERS angstrom\n')
        for r in M:
            f.write('  %.10f  %.10f  %.10f\n' % tuple(r))
        f.write('ATOMIC_POSITIONS crystal\n')
        for e, x, y, z in fr:
            free = '  1 1 1' if e == 'I' else '  0 0 0'
            f.write(f'  {e:<2s}  {x:.10f}  {y:.10f}  {z:.10f}{free}\n')
        f.write('K_POINTS automatic\n  4 4 2 0 0 0\n')
    nfree = sum(1 for e, *_ in fr if e == 'I')
    print(f'{tag}: rerelax_{tag}.in written ({nfree} free iodines, V={abs(np.linalg.det(M)):.1f} A^3)')
