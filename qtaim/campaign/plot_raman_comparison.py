#!/usr/bin/env python3
"""Predicted polyiodide Raman: Y-MOF strained vs activated.
FD phonons on I sublattice (framework frozen), bond-polarizability activities.
Brewer Dark2, solid lines, no title, x-axis = Raman shift, to 400 cm-1.
Tweak freely: COLORS, LW, W (broadening), XLIM, the two .npy paths.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from cycler import cycler

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 8,
    "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.linewidth": 1.0, "lines.linewidth": 1.8,
    "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none",
})

BASE = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
# ---- tweakable parameters ----
XLIM = (40, 225)        # x-axis range
W = 5.0               # Lorentzian broadening half-width (cm-1)
LW = 1.8              # line width
# Brewer Dark2 first two colors
DARK2 = ["#7570B3", "#D95F02"]  # purple, orange
LABELS = ["Y-MOF$_{st}$", "Y-MOF$_{act}$"]
FILES = [
    (f"{BASE}/Y_str/raman/freqs.npy",  f"{BASE}/Y_str/raman/raman_activities.npy"),
    (f"{BASE}/Y_actc/raman/freqs.npy", f"{BASE}/Y_actc/raman/raman_activities.npy"),
]
OUT = f"{BASE}/raman_comparison_Y_str_vs_Y_actc"

def load(freq_file, act_file):
    f = np.load(freq_file); a = np.load(act_file)
    mask = f > 5  # drop imaginary / near-zero
    return f[mask], a[mask]

def broaden(freqs, acts, x, w=W):
    y = np.zeros_like(x)
    for fi, ai in zip(freqs, acts):
        y += ai * w**2 / ((x - fi)**2 + w**2)
    return y

x = np.linspace(XLIM[0], XLIM[1], 1000)
fig, ax = plt.subplots(figsize=(3.33, 2.5))

for i, ((ff, af), color, label) in enumerate(zip(FILES, DARK2, LABELS)):
    freqs, acts = load(ff, af)
    y = broaden(freqs, acts, x)
    y /= y.max()  # normalize
    ax.plot(x, y, color=color, lw=LW, label=label)

ax.set_xlim(*XLIM)
ax.set_xlabel("Raman shift, cm$^{-1}$")
ax.set_ylabel("Raman intensity (norm.)")
ax.legend(loc="upper right")
plt.tight_layout(pad=0.2)
fig.savefig(OUT + ".png", dpi=600)
fig.savefig(OUT + ".pdf")
plt.close()
print(f"saved {OUT}.png/.pdf  (XLIM={XLIM}, W={W}, Dark2 solid)")
