#!/usr/bin/env python3
"""Animations for the headline Raman modes of all five systems.
Same recipe as the approved vc 167.8 test GIF: image-spanning connected motif
(graph walk over the bonded pairs), eigenvector-locked movers, smooth arrows,
distance-based bond drawing over the replicated cells, no framework."""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

C = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442'
OUT = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results_442/raman_442'
AMP = 0.175
NF = 25
DUR = 60

GIFS = [
    ('vc_poly_lS3', 140.2, 'vc_140.2_chain-triangle_stretch'),
    ('vc_poly_lS3', 100.5, 'vc_100.5_secondary-bond_stretch'),
    ('st_poly_lS5', 148.9, 'st_148.9_chain_stretch'),
    ('st_poly_lS5',  89.9, 'st_89.9_link_elongation'),
    ('act_poly_lS2', 158.9, 'act_158.9_branched_unit_stretch'),
    ('act_poly_lS2',  90.4, 'act_90.4_link_stretch'),
    ('act_mol_rd', 211.0, 'actmol_211.0_molecular_stretches'),
    ('adduct', 175.7, 'adduct_175.7_I3_unit_stretch'),
]
BOND_PAIRS = {
    'st_poly_lS5': [(72,73),(73,74),(74,75),(76,77),(77,78),(78,79)],
    'vc_poly_lS3': [(72,76),(72,73),(73,74),(73,75),(76,77),(77,78),(77,79),(78,79)],
    'act_poly_lS2': [(72,73),(74,75),(76,77),(78,79),(73,74),(73,75),(77,78),(77,79)],
    'act_mol_rd': [(72,73),(74,77),(78,79),(75,76)],
    'adduct': [(72,79),(72,73),(77,78),(77,74),(75,76)],
}
NAMES = {'st_poly_lS5': 'Y-MOFst(poly)', 'vc_poly_lS3': 'Y-MOFvc(poly)',
         'act_poly_lS2': 'Y-MOFact(poly)', 'act_mol_rd': 'Y-MOFact(mol)',
         'adduct': 'Y-MOFadduct(test)'}

def make_gif(tag, target, outfile):
    f = np.load(f'{C}/{tag}/raman/freqs.npy')
    evecs = np.load(f'{C}/{tag}/raman/evecs.npy')
    ix = int(np.argmin(np.abs(f - target)))
    e = evecs[:, ix]
    dcart = e / np.sqrt(126.904)
    dcart = dcart / np.abs(dcart).max() * AMP
    diod = dcart.reshape(8, 3)
    movers = {72 + i for i in np.argsort(np.linalg.norm(diod, axis=1))[-2:].tolist()}
    emax = np.linalg.norm(diod, axis=1).max()

    s = open(f'{C}/{tag}/scf_loaded.in').read()
    cl = s.split('CELL_PARAMETERS angstrom\n')[1].split('ATOMIC_POSITIONS')[0].strip().splitlines()
    M = np.array([[float(x) for x in l.split()] for l in cl])
    Minv = np.linalg.inv(M)
    seg = s.split('ATOMIC_POSITIONS crystal\n')[1].split('K_POINTS')[0]
    els, fr = [], []
    for l in seg.strip().splitlines():
        q = l.split()
        els.append(q[0]); fr.append(np.array([float(x) for x in q[1:4]]))
    els = np.array(els); fr = np.array(fr)

    # connected image-spanning placement (BFS over the bonded graph, MIC)
    shifts = {72: np.zeros(3, dtype=int)}
    adj = {}
    for a, b in BOND_PAIRS[tag]:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    queue = [72]
    while queue:
        a = queue.pop(0)
        for b in adj[a]:
            if b in shifts:
                continue
            raw = fr[b] - fr[a]
            n = np.round(raw @ M @ Minv).astype(int)
            shifts[b] = shifts[a] - n
            queue.append(b)

    # disconnected components (separate chains/molecules) -> in-cell placement
    for k in range(72, 80):
        shifts.setdefault(k, np.zeros(3, dtype=int))

    A1 = M[0]; C3 = M[2]
    placed = []
    for n1 in (-1, 0, 1):
        for n3 in (0, 1):
            t = n1 * A1 + n3 * C3
            for k in range(72, 80):
                placed.append(((fr[k] + shifts[k]) @ M + t, k, t))
    P = np.array([p for p, _, _ in placed])
    bond_lines = [(i, j) for i in range(len(placed))
                  for j in range(i + 1, len(placed))
                  if np.linalg.norm(P[i] - P[j]) < 3.45]

    fig = plt.figure(figsize=(7.4, 10.8))
    ax = fig.add_subplot(111, projection='3d')
    frames = []
    for ph in range(NF):
        t = 2*np.pi*ph/NF
        ax.cla()
        icart = np.array([p + diod[k - 72]*np.sin(t) for p, k, _ in placed])
        for i, j in bond_lines:
            ax.plot(*zip(icart[i], icart[j]), color='#2e8b57', lw=1.3, alpha=.85)
        for i, (p, k, _) in enumerate(placed):
            amp_k = np.linalg.norm(diod[k - 72])
            hl = amp_k >= 0.30 * emax and (i // 8 == 1)
            ax.scatter(icart[i,0], icart[i,1], icart[i,2],
                       s=120 if hl else 70, c='#d62728' if hl else '#6a3fb5',
                       edgecolors='k', lw=.4, depthshade=False)
        for i, (p, k, _) in enumerate(placed):
            amp_k = np.linalg.norm(diod[k - 72])
            if amp_k >= 0.30 * emax and (i // 8 == 1):
                d = diod[k - 72]*np.sin(t)
                if np.linalg.norm(d) > 0.008:
                    ax.quiver(icart[i,0], icart[i,1], icart[i,2], d[0]*4, d[1]*4, d[2]*4,
                              color='#e05a2b', lw=2.2, arrow_length_ratio=0.25)
        for i, (p, k, _) in enumerate(placed):
            ax.text(icart[i,0] + 0.32, icart[i,1] + 0.32, icart[i,2], f'I{k - 71}',
                    fontsize=7, color='0.1', ha='center', va='center')
        ax.set_title(f'{NAMES[tag]}: {f[ix]:.1f} cm$^{{-1}}$', fontsize=10)
        ax.set_axis_off()
        ax.view_init(elev=12, azim=-60)
        lo, hi = P.min(0), P.max(0)
        pad = np.array([1.3, 1.3, 1.0])
        ax.set_xlim(lo[0]-pad[0], hi[0]+pad[0])
        ax.set_ylim(lo[1]-pad[1], hi[1]+pad[1])
        ax.set_zlim(lo[2]-pad[2], hi[2]+pad[2])
        fn = f'/tmp/ga_{ph:02d}.png'
        fig.savefig(fn, dpi=110)
        frames.append(Image.open(fn).convert('P', palette=Image.ADAPTIVE))
    frames[0].save(f'{OUT}/{outfile}.gif', save_all=True, append_images=frames[1:],
                   duration=DUR, loop=0)
    print(f'{outfile}.gif: mode {f[ix]:.2f} cm-1, movers ' +
          str(sorted({k - 71 for k in movers})), flush=True)

for tag, target, outfile in GIFS:
    make_gif(tag, target, outfile)
print('ALL GIFS DONE')
