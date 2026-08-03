"""Two-panel publication composite of Multiwfn QTAIM gradient+contour maps.

This is pure image placement: the chemistry (density contours, gradient-line
bond paths, atoms, bond sticks) is rendered entirely by Multiwfn from the cube
files. Here we only trim whitespace, place the two maps side-by-side, and add
panel labels + a caption. No re-rendering of the density field.
"""
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = "Y_act_gradient.png"       # panel (a): covalent I2  (I76-I77, 2.692 A)
B = "Y_deep_gradient.png"      # panel (b): polyiodide   (I75-I74-I76, ~2.93 A)


def trim(im, pad=10, thresh=250):
    """Crop uniform near-white margins, keep a small pad so axis labels survive."""
    g = np.array(im.convert("L"))
    mask = g < thresh                       # non-white pixels
    if not mask.any():
        return im
    ys, xs = np.where(mask)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad, g.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad, g.shape[1])
    return im.crop((x0, y0, x1, y1))


ia = trim(Image.open(A).convert("RGB"))
ib = trim(Image.open(B).convert("RGB"))
# common height
H = min(ia.height, ib.height)
ia = ia.resize((int(ia.width * H / ia.height), H))
ib = ib.resize((int(ib.width * H / ib.height), H))

fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6))
for ax, img, lab, title in [
    (axes[0], ia, "(a)  Y_act — covalent I$_2$ (I76–I77, 2.69 Å)", ""),
    (axes[1], ib, "(b)  Y_deep — polyiodide chain (I75–I74–I76, ~2.93 Å)", ""),
]:
    ax.imshow(np.asarray(img))
    ax.set_axis_off()
    ax.set_title(lab, fontsize=12, pad=6)
fig.suptitle("QTAIM electron density: bond paths and critical points (valence density, ρ)",
             fontsize=13.5, y=0.98)
fig.tight_layout(rect=[0, 0.02, 1, 0.95])
fig.text(0.5, 0.015,
         "Multiwfn gradient-line map from QE cube (val.cube); "
         "solid lines = density contours, trajectories = bond paths, "
         "blue sticks = empirical bonds, axes in Å.",
         ha="center", fontsize=8.5, color="0.3")
fig.savefig("QTAIM_two_panel.png", dpi=300, bbox_inches="tight")
fig.savefig("QTAIM_two_panel.pdf", bbox_inches="tight")
print("wrote QTAIM_two_panel.png / .pdf")
print("panel sizes (px):", ia.size, ib.size)
