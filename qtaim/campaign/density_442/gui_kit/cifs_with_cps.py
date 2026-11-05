#!/usr/bin/env python3
"""Emit kit CIFs with CP dots as He pseudo-atoms.

Usage: python3 cifs_with_cps.py <tag> <cps_txt> <out_cif>
Reads the bare structure CIF Y-MOF<tag>.cif (same dir), appends He sites at the
CP positions from CPs.txt (Bohr Cartesian of the unit cell) and writes <out_cif>.
Run from gui_kit. He dots inherit VESTA's default grey spheres; CP001.. naming.
"""
import sys
import numpy as np

BOHR = 0.5291772109
CELLS = {  # rows a1,a2,a3 in Angstrom (relaxed 442 cells)
    'st':     np.array([[8.9348, 0, 0], [0, 9.6608, 0], [-3.7625602248, 0, 16.5958256665]]),
    'vc':     np.array([[8.7975, 0, 0], [-1.4362, 8.2029, 0], [-4.1207, -0.9336, 17.2148]]),
    'act':    np.array([[8.9348, 0, 0], [0, 9.6608, 0], [-3.7625602248, 0, 16.5958256665]]),
    'actmol': np.array([[8.9348, 0, 0], [0, 9.6608, 0], [-3.7625602248, 0, 16.5958256665]])}

tag, cps_txt, out_cif = sys.argv[1], sys.argv[2], sys.argv[3]
Minv = np.linalg.inv(CELLS[tag] / BOHR)                # bohr Cartesian -> fractional
lines = [l.rstrip('\n') for l in open(f'Y-MOF{tag}.cif')]
out = []
for l in lines:
    out.append(l)
    q = l.split()
    if q and q[0] == 'I8':                             # append after the last atom line
        for cp in open(cps_txt).readlines()[1:]:
            p = cp.split()
            fr = np.array([float(x) for x in p[1:4]]) @ Minv
            out.append(f'CP{int(p[0]):03d} He {fr[0] % 1:.5f} {fr[1] % 1:.5f} {fr[2] % 1:.5f}')
open(out_cif, 'w').write('\n'.join(out) + '\n')
n = sum(1 for l in out if len(l.split()) > 1 and l.split()[1] == 'He')
print(f'{out_cif}: {n} He dots')
