#!/usr/bin/env python3
"""Animation of the vc 167.8 cm-1 bridging-I1-I5 stretch, drawn as a CONNECTED
image-spanning polyiodide motif: each iodine is plotted at the periodic image
that makes its bonded pairs real (graph walk over the 9 bonded pairs, MIC
relative to the already-placed neighbor). The eigenvector displacement of the
image = that of its parent atom. Arrows scale smoothly with sin(t) (no
on/off flicker). Framework dots drawn for the home cell and for each image
cell spanned by the iodines."""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

C = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442/vc_poly_lS3'
OUT = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results_442/raman_442'
TARGET = 167.8
AMP = 0.175
NF = 25
BOND_PAIRS = [(72, 73), (73, 74), (73, 75), (74, 75),      # chain A triangle I1-I2-I3-I4
              (76, 77), (77, 78), (77, 79), (78, 79),      # chain B triangle I5-I6-I7-I8
              (72, 76)]                                     # the bridging I1-I5 molecule

s = open(f'{C}/scf_loaded.in').read()
cl = s.split('CELL_PARAMETERS angstrom\n')[1].split('ATOMIC_POSITIONS')[0].strip().splitlines()
M = np.array([[float(x) for x in l.split()] for l in cl])
Minv = np.linalg.inv(M)
seg = s.split('ATOMIC_POSITIONS crystal\n')[1].split('K_POINTS')[0]
els, fr = [], []
for l in seg.strip().splitlines():
    q = l.split()
    els.append(q[0]); fr.append(np.array([float(x) for x in q[1:4]]))
els = np.array(els); fr = np.array(fr)

f = np.load(f'{C}/raman/freqs.npy')
evecs = np.load(f'{C}/raman/evecs.npy')
ix = int(np.argmin(np.abs(f - TARGET)))
e = evecs[:, ix]
dcart = (e / np.sqrt(126.904))
dcart = dcart / np.abs(dcart).max() * AMP
diod = dcart.reshape(8, 3)                                # iodines only (framework frozen)
_mov = np.argsort(np.linalg.norm(diod, axis=1))[-2:].tolist()
_mov = np.argsort(np.linalg.norm(diod, axis=1))[-2:].tolist()
movers = {72 + i for i in _mov}
print('moving iodines (table labels):', [f'I{i-71}' for i in sorted(movers)])

# ---- image-spanning layout: connected motif over the bonded graph ----
shifts = {72: np.zeros(3, dtype=int)}                       # iodine idx (cube numbering) -> cell shift
adj = {}
for a, b in BOND_PAIRS:
    adj.setdefault(a, []).append(b)
    adj.setdefault(b, []).append(a)
queue = [72]
while queue:
    a = queue.pop(0)
    for b in adj[a]:
        if b in shifts:
            continue
        raw = fr[b] - fr[a]
        n = np.round(raw @ M @ Minv).astype(int)   # nearest lattice vector (frac)
        shifts[b] = shifts[a] - n
        queue.append(b)

# ---- replicate over several cells so the chains are visible ----
A1 = M[0]; C3 = M[2]
TRANSLATIONS = []
for n1 in (-1, 0, 1):
    for n3 in (0, 1):
        TRANSLATIONS.append(n1 * A1 + n3 * C3)
placed = []                                                # (cartesian pos, iodine idx, translation vector)
for t in TRANSLATIONS:
    for k in range(72, 80):
        p = (fr[k] + shifts[k]) @ M + t
        placed.append((p, k, t))
P = np.array([p for p, _, _ in placed])

# ---- distance-based bond detection among placed iodines ----
bond_lines = []
for i in range(len(placed)):
    for j in range(i + 1, len(placed)):
        d = np.linalg.norm(P[i] - P[j])
        if d < 3.45:
            bond_lines.append((i, j))

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
        hl = (k in movers) and (i // 8 == 1)          # highlight movers in the central copy
        ax.scatter(icart[i,0], icart[i,1], icart[i,2],
                   s=120 if hl else 70, c='#d62728' if hl else '#6a3fb5',
                   edgecolors='k', lw=.4, depthshade=False)
    for i, (p, k, _) in enumerate(placed):
        ax.text(icart[i,0] + 0.55, icart[i,1] + 0.55, icart[i,2],
                f'I{k - 71}', fontsize=7.5, color='0.15')
        if (k in movers) and (i // 8 == 1):           # arrows: central copy movers
            d = diod[k - 72]*np.sin(t)
            if np.linalg.norm(d) > 0.01:
                ax.quiver(icart[i,0], icart[i,1], icart[i,2], d[0]*4, d[1]*4, d[2]*4,
                          color='#e05a2b', lw=2.2, arrow_length_ratio=0.25)
    ax.set_title(f'Y-MOFvc(poly): {f[ix]:.1f} cm$^{{-1}}$ — bridging I1–I5 stretch', fontsize=10)
    ax.set_axis_off()
    ax.view_init(elev=12, azim=-60)
    fn = f'/tmp/gif_{ph:02d}.png'
    fig.savefig(fn, dpi=110)
    frames.append(Image.open(fn).convert('P', palette=Image.ADAPTIVE))
frames[0].save(f'{OUT}/vc_167.8_bridged_I2_stretch.gif', save_all=True,
               append_images=frames[1:], duration=60, loop=0)
print(f'GIF written ({NF} frames, 6 copies, {len(bond_lines)} bonds)')
