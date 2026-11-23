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
AMP = 0.35
NF = 13
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

# ---- image-spanning layout: walk the bonded graph, MIC each new atom ----
shifts = {72: np.zeros(3, dtype=int)}                       # iodine idx -> cell shift
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
disp_pos = {k: (fr[k] + shifts.get(k, np.zeros(3, dtype=int))) @ M for k in range(72, 80)}
d_eig = {k: diod[k - 72] for k in range(72, 80)}                     # image = parent atom's eigenvector

# framework dots: home cell + each spanned image cell
fw_pts = []
fw_idx = [k for k in range(len(els)) if els[k] != 'I']
uniq_shifts = {tuple(v) for v in shifts.values()}
for sh in uniq_shifts:
    sh = np.array(sh)
    for k in fw_idx:
        fw_pts.append((fr[k] + sh) @ M)
fw_pts = np.array(fw_pts)

fig = plt.figure(figsize=(7.2, 6.4))
ax = fig.add_subplot(111, projection='3d')
frames = []
for ph in range(NF):
    t = 2*np.pi*ph/NF
    ax.cla()
    ax.scatter(fw_pts[:,0], fw_pts[:,1], fw_pts[:,2], s=14, c='0.55', depthshade=False)
    icart = {k: disp_pos[k] + d_eig[k]*np.sin(t) for k in disp_pos}
    for k in disp_pos:
        hl = k in movers and abs(np.sin(t)) > 0.05
        ax.scatter(icart[k][0], icart[k][1], icart[k][2],
                   s=210 if hl else 150, c='#d62728' if hl else '#6a3fb5',
                   edgecolors='k', lw=.6, depthshade=False)
    for k in disp_pos:
        d = d_eig[k]*np.sin(t)
        n = np.linalg.norm(d)
        if n > 0.005:
            ax.quiver(icart[k][0], icart[k][1], icart[k][2], d[0]*6, d[1]*6, d[2]*6,
                      color='#e05a2b', lw=2.4, arrow_length_ratio=0.25)
    for a, b in BOND_PAIRS:
        ax.plot(*zip(icart[a], icart[b]), color='#2e8b57', lw=1.4, alpha=.85)
    ax.set_title(f'Y-MOFvc(poly): {f[ix]:.1f} cm$^{{-1}}$ — bridging I1–I5 stretch\n'
                 f'red = moving I1/I5; arrows = displacements (exaggerated); '
                 f'green = bonded I-I pairs; the motif spans a cell edge', fontsize=9.5)
    ax.set_axis_off()
    ax.view_init(elev=12, azim=-60)
    fn = f'/tmp/gif_{ph:02d}.png'
    fig.savefig(fn, dpi=110)
    frames.append(Image.open(fn).convert('P', palette=Image.ADAPTIVE))
frames[0].save(f'{OUT}/vc_167.8_bridged_I2_stretch.gif', save_all=True,
               append_images=frames[1:], duration=90, loop=0)
print(f'GIF written ({NF} frames); shifts: ' +
      ', '.join(f'I{k+1}:{tuple(v)}' for k, v in shifts.items()))
