#!/usr/bin/env python3
"""Animation test (user request 3): GIF of the vc 167.8 cm-1 bridging-I1-I5
stretch. Iodine displacement arrows (Chemcraft-style) + oscillating iodine
positions against the frozen framework. PIL-assembled animated GIF."""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
sys.path.insert(0, '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/bcp_plots/multiwfn_out')
from seed_bcp import load_cube, BOHR

C = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442/vc_poly_lS3'
OUT = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results_442/raman_442'
TAGET = 167.8
AMP = 0.35          # max displacement amplitude (A)
NF = 13

f = np.load(f'{C}/raman/freqs.npy')
evecs = np.load(f'{C}/{ "raman" }/evecs.npy')
ix = int(np.argmin(np.abs(f - TAGET)))
e = evecs[:, ix]
dcart = (e / np.sqrt(126.904))            # mass-weighted -> Cartesian (arb units)
dcart = dcart / np.abs(dcart).max() * AMP # normalize max atomic amplitude

s = open(f'{C}/scf_loaded.in').read()
cl = s.split('CELL_PARAMETERS angstrom\n')[1].split('ATOMIC_POSITIONS')[0].strip().splitlines()
M = np.array([[float(x) for x in l.split()] for l in cl])
seg = s.split('ATOMIC_POSITIONS crystal\n')[1].split('K_POINTS')[0]
els, fr = [], []
for l in seg.strip().splitlines():
    q = l.split()
    els.append(q[0]); fr.append(np.array([float(x) for x in q[1:4]]))
els = np.array(els); fr = np.array(fr)
Minv = np.linalg.inv(M)
cart = fr @ M
dc = dcart.reshape(8, 3) if dcart.size == 24 else dcart[:24].reshape(8, 3)
# pad if needed (framework atoms have zero displacement anyway)
if dcart.size == 24:
    dfull = np.zeros((80, 3)); dfull[72:] = dc
else:
    dfull = dcart

fig = plt.figure(figsize=(7.2, 6.4))
ax = fig.add_subplot(111, projection='3d')
# movers fixed ONCE from the eigenvector (phase-independent): at the
# equilibrium-crossing frames the per-frame amplitudes vanish into noise and
# the highlight would flicker onto the symmetry-equivalent chain iodines.
movers = set(np.argsort(np.linalg.norm(dcart.reshape(8, 3), axis=1))[-2:].tolist())
frames = []
for ph in range(NF):
    t = 2*np.pi*ph/NF
    disp = dfull * np.sin(t)
    ax.cla()
    fw = els != 'I'
    ax.scatter(cart[fw,0], cart[fw,1], cart[fw,2], s=14, c='0.55', depthshade=False)
    icart = cart[72:] + disp[72:]
    for k in range(8):
        hl = k in movers
        ax.scatter(icart[k,0], icart[k,1], icart[k,2],
                   s=200 if hl else 150, c='#d62728' if hl else '#6a3fb5',
                   edgecolors='k', lw=.6, depthshade=False)
    for k in range(8):
        d = disp[72+k]
        n = np.linalg.norm(d)
        if n > 0.03:
            ax.quiver(icart[k,0], icart[k,1], icart[k,2], d[0]*6, d[1]*6, d[2]*6,
                      color='#e05a2b', lw=2.4, arrow_length_ratio=0.25)
    for i, j in [(72,73),(73,74),(73,75),(74,75),(76,77),(77,78),(77,79),(78,79),(72,76)]:
        ax.plot(*zip(cart[i], cart[j]), color='#2e8b57', lw=1.4, alpha=.85)
    ax.set_title(f'Y-MOFvc(poly): {f[ix]:.1f} cm$^{{-1}}$ — bridging I1–I5 stretch\n'
                 f'red = the two moving iodines; arrows = displacement vectors (exaggerated); '
                 f'green = bonded I-I pairs', fontsize=10)
    ax.set_axis_off()
    ax.view_init(elev=12, azim=-60)
    fn = f'/tmp/gif_{ph:02d}.png'
    fig.savefig(fn, dpi=110)
    frames.append(Image.open(fn).convert('P', palette=Image.ADAPTIVE))
frames[0].save(f'{OUT}/vc_167.8_bridged_I2_stretch.gif', save_all=True,
               append_images=frames[1:], duration=90, loop=0)
print(f'GIF written: vc_167.8_bridged_I2_stretch.gif ({NF} frames, mode {f[ix]:.1f} cm-1)')
