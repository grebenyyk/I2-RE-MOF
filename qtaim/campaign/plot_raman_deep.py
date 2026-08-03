#!/usr/bin/env python3
"""Deep-minimum Raman (sandwiched I2 + uniform chain): Gd/Tb/Y overlaid.
FD phonons on I sublattice (framework frozen), bond-polarizability activities.
Shows the metal-independent signature: ~225 (I2 stretch) + ~241-247 (chain multiplet).
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 8,
    "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.linewidth": 1.0, "lines.linewidth": 1.8,
    "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none",
})

BASE = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
XLIM = (40, 225)
W = 5.0           # Lorentzian broadening half-width (cm-1)
LW = 1.8
# 3 distinct, colorblind-friendly colors
COL = ["#1b9e77", "#d95f02", "#7570b3"]   # teal, orange, purple
LABELS = ["Gd-MOF$_{\\rm deep}$", "Tb-MOF$_{\\rm deep}$", "Y-MOF$_{\\rm deep}$"]
FILES = [
    (f"{BASE}/Gd_deep/raman/freqs.npy", f"{BASE}/Gd_deep/raman/raman_activities.npy"),
    (f"{BASE}/Tb_deep/raman/freqs.npy", f"{BASE}/Tb_deep/raman/raman_activities.npy"),
    (f"{BASE}/Y_deep/raman/freqs.npy",  f"{BASE}/Y_deep/raman/raman_activities.npy"),
]
OUT = f"{BASE}/raman_deep_Gd_Tb_Y"

def load(ff, af):
    f = np.load(ff); a = np.load(af)
    m = f > 5
    return f[m], a[m]

def broaden(freqs, acts, x, w=W):
    y = np.zeros_like(x)
    for fi, ai in zip(freqs, acts):
        y += ai * w**2 / ((x - fi)**2 + w**2)
    return y

x = np.linspace(XLIM[0], XLIM[1], 1200)
fig, ax = plt.subplots(figsize=(3.33, 2.5))

for (ff, af), color, label in zip(FILES, COL, LABELS):
    freqs, acts = load(ff, af)
    y = broaden(freqs, acts, x)
    y /= y.max()
    ax.plot(x, y, color=color, lw=LW, label=label)

ax.set_xlim(*XLIM)
ax.set_xlabel("Raman shift (cm$^{-1}$)")
ax.set_ylabel("Raman intensity (norm.)")
# annotate the I2 stretch (~225) and chain multiplet (~245)
ax.axvline(213, color="grey", ls=":", lw=0.8)
ax.text(213, 1.02, "gas I$_2$", color="grey", fontsize=6.5, ha="center", transform=ax.get_xaxis_transform())
ax.legend(loc="upper left")
plt.tight_layout(pad=0.2)
fig.savefig(OUT + ".png", dpi=600)
fig.savefig(OUT + ".pdf")
plt.close()
print(f"saved {OUT}.png/.pdf  (XLIM={XLIM}, W={W}, 3 deep metals overlaid)")
