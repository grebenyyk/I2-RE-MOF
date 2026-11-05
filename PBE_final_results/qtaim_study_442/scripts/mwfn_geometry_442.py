#!/usr/bin/env python3
"""Geometry groundwork for the 442 Multiwfn masters (MIC-corrected).

Cube atom order reproduces the table labels: atoms 73..80 = I1..I8 (all
BCP_II_full_442 d_pair values reproduced). Here: for candidate 3-atom planes,
the MIC-corrected perpendicular distance of every iodine to the plane, scored
by how many table-listed I-I bonds have BOTH endpoints in the plane (|d|<0.8 A)
— the bonds whose bond paths will actually draw in the 2D slice.
"""
import sys
import numpy as np
sys.path.insert(0, '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/bcp_plots/multiwfn_out')
from seed_bcp import load_cube, BOHR

CAMP = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442'
SYS = [('st', 'st_vc', 'st_poly_lS5'), ('vc', 'st_vc', 'vc_poly_lS3'),
       ('act', 'act', 'act_poly_lS2'), ('actmol', 'act', 'act_mol_rd')]

# table-listed I-I bonds (labels) per system, for the in-plane coverage score
BONDS = {
    'st': [(1, 2), (2, 3), (3, 4), (5, 6), (6, 7), (7, 8), (1, 4), (5, 8), (1, 3), (5, 7)],
    'vc': [(1, 5), (2, 3), (3, 4), (2, 4), (6, 7), (6, 8), (7, 8), (1, 2), (5, 6)],
    'act': [(1, 2), (2, 3), (2, 4), (5, 6), (6, 8), (6, 7), (3, 4), (7, 8)],
    'actmol': [(1, 2), (3, 6), (4, 5), (7, 8)],
}
CAND = {
    'st': [(1, 2, 3), (1, 2, 4), (5, 6, 7)],
    'vc': [(1, 2, 5), (1, 2, 6), (1, 5, 3), (1, 5, 6), (2, 3, 4)],
    'act': [(1, 2, 3), (1, 2, 4)],
    'actmol': [(1, 2, 3), (1, 2, 5), (1, 2, 6), (4, 5, 7)],
}

for tag, d, name in SYS:
    origin, H, dims, Z, pos, rho = load_cube(f'{CAMP}/{d}/rho17_{name}.cube')
    M = H.T * np.array(dims, dtype=float)[:, None]
    Minv = np.linalg.inv(M)
    ilab = {k + 1: pos[72 + k] for k in range(8)}     # label -> Bohr Cartesian

    def micv(a, b):
        """MIC displacement vector a->b (Bohr)."""
        dd = (ilab[b] - ilab[a]) @ Minv
        return (dd - np.round(dd)) @ M

    print(f'\n=== {tag} ===')
    best = None
    for tri in CAND[tag]:
        P = ilab[tri[0]]
        v1, v2 = micv(tri[0], tri[1]), micv(tri[0], tri[2])
        n = np.cross(v1, v2); n /= np.linalg.norm(n)
        d = {L: abs(np.dot(micv(tri[0], L), n)) * BOHR for L in ilab}
        cov = [(a, b) for (a, b) in BONDS[tag] if d[a] < 0.8 and d[b] < 0.8]
        missed = [(a, b) for (a, b) in BONDS[tag] if not (d[a] < 0.8 and d[b] < 0.8)]
        trilab = ','.join(f'I{t}' for t in tri)
        print(f'  plane {trilab}:  in-plane bonds {len(cov)}/{len(BONDS[tag])} -> {cov}')
        print(f'     out-of-plane: ' + '  '.join(f'I{L}:{d[L]:.2f}' for L in sorted(d) if L not in tri))
        if missed:
            print(f'     misses: {missed}')
        if best is None or len(cov) > best[0]:
            best = (len(cov), tri, d)
    _, tri, d = best
    print(f'  -> BEST: plane I{tri[0]}-I{tri[1]}-I{tri[2]}')
