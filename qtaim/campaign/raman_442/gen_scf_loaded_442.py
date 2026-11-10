#!/usr/bin/env python3
"""After the tight I-only re-relax converges, build <sys>/scf_loaded.in: an SCF
at the FINAL re-relaxed geometry, prefix raman442_<sys>_loaded, outdir ./tmp_loaded
(the raman_driver_442 seeding expects exactly this). The converged re-relax .save
is copied into tmp_loaded so the SCF restarts from the converged wavefunction
(same geometry, no basin trap). Usage: gen_scf_loaded_442.py <tag>"""
import re, shutil, sys, os

BASE = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442'
tag = sys.argv[1]
d = f'{BASE}/{tag}'
out = f'{d}/rerelax_{tag}.out'
txt = open(out, errors='ignore').read()
assert 'JOB DONE' in txt, f'{tag}: re-relax not finished'
tail = txt[txt.index('Begin final coordinates'):]
cm = re.search(r'CELL_PARAMETERS \(angstrom\)\s*\n\s*(.+)\n\s*(.+)\n\s*(.+)\n', tail)
cell = [m.split() for m in cm.groups()]
pm = re.search(r'ATOMIC_POSITIONS \(crystal\)\s*\n((?:\s*\S+\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+\s+[01]\s+[01]\s+[01].*\n)+)', tail)
pos = [l.split() for l in pm.group(1).splitlines()]
assert len(pos) == 80

src_save = f'{d}/tmp_rerelax/raman442_{tag}_rerelax.save'
assert os.path.isdir(src_save), f'missing re-relax save: {src_save}'
os.makedirs(f'{d}/tmp_loaded', exist_ok=True)
dst_save = f'{d}/tmp_loaded/raman442_{tag}_loaded.save'
if not os.path.isdir(dst_save):
    shutil.copytree(src_save, dst_save)

head = open(f'{d}/rerelax_{tag}.in').read().split('ATOMIC_SPECIES')[0]
head = head.replace("calculation = 'relax'", "calculation = 'scf'")
head = head.replace("prefix = 'raman442_%s_rerelax'" % tag, "prefix = 'raman442_%s_loaded'" % tag)
head = head.replace("outdir = './tmp_rerelax'", "outdir = './tmp_loaded'")
head = head.replace('  etot_conv_thr = 1.0d-6\n  forc_conv_thr = 1.0d-5\n  nstep = 300\n', '  nstep = 1\n')
head = head.replace('&IONS\n  ion_dynamics = \'bfgs\'\n/', '')
if 'restart_mode' not in head:
    head = head.replace("  calculation = 'scf'", "  calculation = 'scf'\n  restart_mode = 'restart'")
species = 'ATOMIC_SPECIES' + open(f'{d}/rerelax_{tag}.in').read().split('ATOMIC_SPECIES')[1].split('CELL_PARAMETERS')[0]

with open(f'{d}/scf_loaded.in', 'w') as f:
    f.write(head)
    f.write(species)
    f.write('CELL_PARAMETERS angstrom\n')
    for r in cell:
        f.write('  ' + '  '.join(r) + '\n')
    f.write('ATOMIC_POSITIONS crystal\n')
    for q in pos:
        f.write('  %s  %s  %s  %s\n' % (q[0], q[1], q[2], q[3]))
    f.write('K_POINTS automatic\n  4 4 2 0 0 0\n')
print(f'{tag}: scf_loaded.in written (final geometry, restart from re-relax save)')
