#!/usr/bin/env python3
"""Final consistency-checked package for the 442 Raman campaign.

Creates PBE_final_results_442/raman_442/structures_and_analysis/ with:
  - five CIFs = the exact geometries the FD spectra were computed on
    (iodine labels I1..I8 = Table S11 / cube order; framework labels = l442 CIF order)
  - energy_decomposition_442.csv (frozen-fragment decomposition, 3 variants)
  - vibration_modes_442.csv (per-mode: band, rel. intensity, iodine atoms,
    natural-language description; assignments from eigenvector projections)

Naming convention (binding, matches Table S11):
  Y-MOFst(poly) / Y-MOFvc(poly) / Y-MOFact(poly) / Y-MOFact(mol) / Y-MOFadduct(test)
  iodines I1..I8 = last 8 atoms of the 80-atom cell, in cell order.
"""
import os, csv, sys
import numpy as np

C = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442'
CIFS = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results_442/landscape_3variants/cifs'
OUT = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results_442/raman_442/structures_and_analysis'
os.makedirs(OUT, exist_ok=True)
BOHR = 0.5291772109

SYS = [('st_poly_lS5',    'Y-MOFst(poly)',        'Y-MOFst_poly.cif',      'l442_st_poly_lS5_relaxed.cif'),
       ('vc_poly_lS3',    'Y-MOFvc(poly)',        'Y-MOFvc_poly.cif',      'l442_vc_poly_lS3_relaxed.cif'),
       ('act_poly_lS2',   'Y-MOFact(poly)',       'Y-MOFact_poly.cif',     'l442_act_poly_lS2_relaxed.cif'),
       ('act_mol_rd',     'Y-MOFact(mol)',        'Y-MOFact_mol.cif',      'l442_act_mol_rd_relaxed.cif'),
       ('adduct',         'Y-MOFadduct(test)',    'Y-MOFadduct_test.cif',  'l442_act_mol_rd_relaxed.cif')]

# ---------- CIFs from the exact FD geometries ----------
def labels_of(cif):
    lab = []
    for l in open(cif):
        q = l.split()
        if len(q) >= 5 and q[0][0].isalpha() and (q[0][1:].isdigit() or q[0] in ('I','C','O','Y','H','N')) and q[1] in ('C','H','N','O','Y','I'):
            if q[0] not in ('_atom_site',) and not q[0].startswith('_'):
                lab.append(q[0])
    return lab

def geom_of(scf):
    s = open(scf).read()
    cl = s.split('CELL_PARAMETERS angstrom\n')[1].split('ATOMIC_POSITIONS')[0].strip().splitlines()
    M = np.array([[float(x) for x in l.split()] for l in cl])
    seg = s.split('ATOMIC_POSITIONS crystal\n')[1].split('K_POINTS')[0]
    fr = []
    for l in seg.strip().splitlines():
        q = l.split()
        if len(q) >= 4 and q[0] in ('C','H','N','O','Y','I'):
            fr.append((q[0], np.array([float(x) for x in q[1:4]])))
    return M, fr

a_, b_, c_ = None, None, None
for tag, name, cifout, labcif in SYS:
    M, fr = geom_of(f'{C}/{tag}/scf_loaded.in')
    labs = labels_of(f'{CIFS}/{labcif}')
    assert len(labs) == len(fr) == 80, (tag, len(labs), len(fr))
    a = np.linalg.norm(M[0]); b = np.linalg.norm(M[1]); c = np.linalg.norm(M[2])
    al = np.degrees(np.arccos(np.dot(M[1],M[2])/(b*c)))
    be = np.degrees(np.arccos(np.dot(M[0],M[2])/(a*c)))
    ga = np.degrees(np.arccos(np.dot(M[0],M[1])/(a*b)))
    with open(f'{OUT}/{cifout}', 'w') as f:
        f.write(f'data_{name}\n'
                f'_audit_creation_method  "QE 7.5 PBE+D3(BJ), 4x4x2-relaxed; iodine sublattice '
                f're-optimized at the Raman production footing; framework frozen"\n')
        f.write(f'_cell_length_a {a:.6f}\n_cell_length_b {b:.6f}\n_cell_length_c {c:.6f}\n')
        f.write(f'_cell_angle_alpha {al:.4f}\n_cell_angle_beta {be:.4f}\n_cell_angle_gamma {ga:.4f}\n')
        f.write("_symmetry_space_group_name_H-M 'P 1'\n_symmetry_Int_Tables_number 1\n")
        f.write('loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n')
        for i, (e, frv) in enumerate(fr):
            f.write(f'{labs[i]} {e} {frv[0]%1:.6f} {frv[1]%1:.6f} {frv[2]%1:.6f}\n')
    print(f'{cifout}: 80 atoms, labels {labs[72:]}')
print('CIFs done')
