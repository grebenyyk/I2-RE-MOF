#!/usr/bin/env python3
"""Recompute Abramov G/V/H and the 4 correlation energies for ALL rows of ALL 442 CP
tables from the final rho_BCP [e/bohr^3] and lap [e/A^5] values, so every row is
self-consistent regardless of which footing produced it. Adds the E columns to the
I-O / I-C / Y-I tables (the I-I table already has them).

Units in the tables: rho(e/A3) = rho_au/BOHR^3; lap(e/A5) = lap_au/BOHR^5;
G,V,H in Hartree/bohr^3 (a.u.); energies in kJ/mol (Hartree x 2625.4996).
Abramov: G = (3/10)(3pi^2)^(2/3) rho^(5/3) + lap/6   [a.u.]
         V = lap/4 - 2G                              [a.u.]
         H = G + V
Correlations: E1 = -V/2 (Espinosa-Molins-Lecomte); E2 = 0.429 G (Abrov/G-based);
              E3 = 0.68 |V| (iodine-specific); E4 = 0.67 G (iodine-specific).
"""
import csv, os
import numpy as np
BOHR = 0.5291772109
BOHR5 = BOHR**5
HA = 2625.4996
TF = 0.3 * (3*np.pi**2)**(2/3)
T = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/qtaim_study_442/tables'
ECOLS = ['E=-V/2', 'E=0.429G', 'E=0.68|V|', 'E=0.67G']
for fname in ('BCP_II_full_442.csv', 'BCP_IO_full_442.csv', 'BCP_IC_full_442.csv', 'BCP_MI_native_442.csv'):
    p = f'{T}/{fname}'
    rows = list(csv.DictReader(open(p)))
    hdr = list(rows[0].keys()) if rows else []
    for col in ECOLS:
        if col not in hdr:
            # insert after H(au) (or before src)
            i = hdr.index('H(au)') + 1 if 'H(au)' in hdr else len(hdr) - 1
            hdr.insert(i, col)
    for r in rows:
        rho = float(r['rho_BCP(e/bohr3)'])
        lap_au = float(r['lap(e/A5)']) * BOHR5
        g = TF * rho**(5/3) + lap_au / 6.0
        v = lap_au / 4.0 - 2.0 * g
        h = g + v
        r['G(au)'] = f'{g:.4f}'
        r['V(au)'] = f'{v:.4f}'
        r['H(au)'] = f'{h:.4f}'
        r['E=-V/2'] = f'{-v/2*HA:.1f}'
        r['E=0.429G'] = f'{0.429*g*HA:.1f}'
        r['E=0.68|V|'] = f'{0.68*abs(v)*HA:.1f}'
        r['E=0.67G'] = f'{0.67*g*HA:.1f}'
    with open(p, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=hdr)
        w.writeheader()
        w.writerows(rows)
    print(f'{fname}: {len(rows)} rows, G/V/H + 4 energies recomputed')
