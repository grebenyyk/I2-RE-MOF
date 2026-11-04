#!/usr/bin/env python3
"""Individual Bader delta-rho charges on the O and I atoms taking part in I-O
non-covalent interactions (pairs with a (3,-1) BCP) vs the I-O distance.
Left pane: charge on O; right pane: charge on I. Systems: st/vc/act/act(mol).
Pairs re-derived with the same logic as the generator (fixed-cell classification:
nearest-2 nuclei, d1<2.6 A, d2/d1<2, ratio-guarded iodine 4d-shell filter, pair={O,I}).
Per-atom charges: ct_{tag}.out 4th column (old verified convention), atom = index in
cube order (critic2 attractor order).
"""
import sys, re
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/bcp_plots/multiwfn_out')
from seed_bcp import load_cube, BOHR

CAMP = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442'
OUT = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/qtaim_study_442'
ROWRE = re.compile(r"\(\s*3\s*,\s*-1\s*\)\s+bond\s+([0-9.+-]+)\s+([0-9.+-]+)\s+([0-9.+-]+)\s+\d+\s+\S+\s+([0-9.Ee+-]+)\s+[0-9.Ee+-]+\s+([0-9.Ee+-]+)")
CTRE = re.compile(r'^\s*(\d+)\s+\d+\s+\d+\s+([A-Z][a-z]?)\s+\d+\s+--\s+([0-9.+\-E]+)\s+([0-9.+\-E]+)\s+([0-9.+\-E]+)\s+([0-9.+\-E]+)')
SYS = [('Y-MOFst', 'st_vc', 'st_poly_lS5'), ('Y-MOFvc', 'st_vc', 'vc_poly_lS3'),
       ('Y-MOFact', 'act', 'act_poly_lS2'), ('Y-MOFact(mol.)', 'act', 'act_mol_rd')]
COLOR = {'Y-MOFst': '#1f77b4', 'Y-MOFvc': '#d62728', 'Y-MOFact': '#2ca02c', 'Y-MOFact(mol.)': '#7f7f7f'}
MARK = {'Y-MOFst': 'o', 'Y-MOFvc': 's', 'Y-MOFact': 'D', 'Y-MOFact(mol.)': '^'}

pairs_all = []
for lbl, d, tag in SYS:
    D = f'{CAMP}/{d}'
    origin, H, dims, Z, pos, rho = load_cube(f'{D}/rho17_{tag}.cube')
    M = H.T * np.array(dims, dtype=float)[:, None]
    Minv = np.linalg.inv(M)
    def mic(a, b):
        dd = (pos[a]-pos[b]) @ Minv; dd -= np.round(dd)
        return np.linalg.norm(dd @ M) * BOHR
    def dp(x):
        return sorted((np.linalg.norm((pos[k]-x-np.round((pos[k]-x)@Minv)@M))*BOHR, k) for k in range(len(Z)))
    ct = {}
    for l in open(f'{D}/ct_{tag}.out', errors='ignore'):
        m = CTRE.match(l)
        if m: ct[int(m.group(1))-1] = float(m.group(6))
    pairs = []
    for l in open(f'{D}/bcp_raw_{tag}.out', errors='ignore'):
        m = ROWRE.search(l)
        if not m: continue
        fr = np.array([float(m.group(i)) for i in (1, 2, 3)])
        x = fr @ M
        ds = dp(x)
        d1, d2, k1, k2 = ds[0][0], ds[1][0], ds[0][1], ds[1][1]
        if not (d1 < 2.6 and d2/d1 < 2): continue
        if (Z[k1] == 53 and d1 < 1.35 and d2/d1 > 1.5) or (Z[k2] == 53 and d2 < 1.35 and d2/d1 > 1.5): continue
        if {Z[k1], Z[k2]} != {53, 8}: continue
        oi = k1 if Z[k1] == 8 else k2
        ii = k1 if Z[k1] == 53 else k2
        pairs.append((mic(oi, ii), ct.get(oi), ct.get(ii)))
    print(f'{lbl}: {len(pairs)} I-O pairs')
    pairs_all.append((lbl, pairs))

fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), sharey=False)
for lbl, pairs in pairs_all:
    if not pairs: continue
    pts = sorted(pairs)
    axes[0].scatter([p[0] for p in pts], [p[1] for p in pts], s=30, color=COLOR[lbl],
                    marker=MARK[lbl], alpha=.85, edgecolors='k', linewidths=.4, label=lbl)
    axes[1].scatter([p[0] for p in pts], [p[2] for p in pts], s=30, color=COLOR[lbl],
                    marker=MARK[lbl], alpha=.85, edgecolors='k', linewidths=.4)
for ax, yl, ttl in ((axes[0], r'$\delta\rho(\mathrm{O})$ / e', 'charge on O'),
                    (axes[1], r'$\delta\rho(\mathrm{I})$ / e', 'charge on I')):
    ax.axhline(0, color='k', lw=.6)
    ax.set_xlabel(r'$d(\mathrm{I}\cdots\mathrm{O})$ / $\mathrm{\AA}$')
    ax.set_ylabel(yl); ax.set_title(ttl, fontsize=10); ax.tick_params(labelsize=8)
axes[0].legend(fontsize=7, frameon=False, loc='lower right')
fig.tight_layout()
fig.savefig(f'{OUT}/charges_vs_IO_distance.png', dpi=200)
fig.savefig(f'{OUT}/charges_vs_IO_distance.pdf')
print('saved charges_vs_IO_distance.png/.pdf')
