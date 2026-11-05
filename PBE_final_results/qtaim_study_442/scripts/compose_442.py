#!/usr/bin/env python3
"""Two-by-two composite of the 442 Multiwfn QTAIM masters with panel letters.

Pure image placement (like the Gamma-era compose_two_panel.py): trims the
rendered maps, places them on a 2x2 grid, adds (a)-(d), titles and a caption.
The chemistry is rendered entirely by Multiwfn from the valence cubes.
Run from qtaim/campaign/density_442/mwfn_out; outputs land in the results dir.
"""
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORK = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/mwfn_out'
OUT = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/qtaim_study_442/multiwfn_out'

PANELS = [
    ('mwfn_442_st.png',     '(a)  Y-MOF$_\\mathregular{st}$ — zigzag polyiodide chain'),
    ('mwfn_442_vc.png',     '(b)  Y-MOF$_\\mathregular{vc}$ — I$_2$ bridging two zigzag chains'),
    ('mwfn_442_act.png',    '(c)  Y-MOF$_\\mathregular{act}$ — branched polyiodide chain'),
    ('mwfn_442_actmol.png', '(d)  Y-MOF$_\\mathregular{act}$(mol.) — molecular I$_2$'),
]


def trim(im, pad=12, thresh=250):
    g = np.array(im.convert('L'))
    mask = g < thresh
    if not mask.any():
        return im
    ys, xs = np.where(mask)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad, g.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad, g.shape[1])
    return im.crop((x0, y0, x1, y1))


ims = [trim(Image.open(f'{WORK}/{f}').convert('RGB')) for f, _ in PANELS]
H = min(im.height for im in ims)
ims = [im.resize((int(im.width * H / im.height), H)) for im in ims]

fig, axes = plt.subplots(2, 2, figsize=(12.6, 11.0))
for ax, (im, (_, title)) in zip(axes.flat, zip(ims, PANELS)):
    ax.imshow(np.asarray(im))
    ax.set_axis_off()
    ax.set_title(title, fontsize=12.5, pad=6)
fig.suptitle('QTAIM electron density of I$_2$@Y-MOF at 4×4×2 $k$-mesh: bond paths and critical points',
             fontsize=14, y=0.985)
fig.text(0.5, 0.012,
         'Multiwfn gradient-line maps from the PAW-valence cubes (val, Γ-cell geometry relaxed at 4×4×2); '
         'solid lines = density contours, trajectories = bond paths, blue dots = (3,−1) critical points '
         '(iodine-involving, ρ ≥ 0.003 a.u.), blue sticks = empirical bonds, red labels = atoms within 0.26 Å '
         'of the plane, axes in Å. Dots follow the validated 442 CP tables; black boundaries mark the periodic-box edge.',
         ha='center', fontsize=8.5, color='0.3')
fig.tight_layout(rect=[0, 0.025, 1, 0.965])
fig.savefig(f'{OUT}/QTAIM_442_two_by_two.png', dpi=300, bbox_inches='tight')
fig.savefig(f'{OUT}/QTAIM_442_two_by_two.pdf', bbox_inches='tight')
print('wrote QTAIM_442_two_by_two.png/.pdf')
print('panel sizes (px):', [im.size for im in ims])
