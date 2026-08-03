#!/usr/bin/env python3
"""Y-MOF Raman, three motifs overlaid (st/act/vc) + gas-phase I2 reference.
FD phonons on I sublattice, bond-polarizability. Frequencies use the corrected
Ry->cm^-1 factor (3634.77); gas I2 validates at 214.5 vs exp 213."""
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
XLIM = (40, 300)
W = 5.0
LW = 1.8
COL = ["#1b9e77", "#d95f02", "#e6ab02"]   # Dark2: teal=st, orange=act, gold=vc (purple reserved for I2)
LABELS = ["Y-MOF$_{\\rm st}$", "Y-MOF$_{\\rm act}$", "Y-MOF$_{\\rm vc}$"]
FILES = [
    (f"{BASE}/Y_str/raman/freqs.npy",  f"{BASE}/Y_str/raman/raman_activities.npy"),
    (f"{BASE}/Y_actc/raman/freqs.npy", f"{BASE}/Y_actc/raman/raman_activities.npy"),
    (f"{BASE}/Y_deep/raman/freqs.npy",  f"{BASE}/Y_deep/raman/raman_activities.npy"),
]
I2_COL = "#7570b3"   # Dark2 purple = gas I2
OUT = f"{BASE}/raman_Y_st_act_vc"

def load(ff, af):
    f = np.load(ff); a = np.load(af)
    m = f > 5
    return f[m], a[m]

def broaden(freqs, acts, x, w=W):
    y = np.zeros_like(x)
    for fi, ai in zip(freqs, acts):
        y += ai * w**2 / ((x - fi)**2 + w**2)
    return y

x = np.linspace(XLIM[0], XLIM[1], 1400)
fig, ax = plt.subplots(figsize=(3.33, 2.5))

for (ff, af), color, label in zip(FILES, COL, LABELS):
    freqs, acts = load(ff, af)
    y = broaden(freqs, acts, x)
    y /= y.max()
    ax.plot(x, y, color=color, lw=LW, label=label)

# gas-phase I2 reference spectrum (validated: stretch 214.5 cm-1 vs exp 213)
fgas = float(np.load("/workspace/dimitry/calc/qe/i2_mof/i2_gas_raman/freqs_gas.npy").max())
ygas = W**2 / ((x - fgas)**2 + W**2); ygas /= ygas.max()
ax.plot(x, ygas, color=I2_COL, lw=LW, label="gas I$_2$")

ax.set_xlim(*XLIM)
ax.set_xlabel("Raman shift, cm$^{-1}$")
ax.set_ylabel("Raman intensity, norm.")
ax.legend(loc="upper left")
plt.tight_layout(pad=0.2)
fig.savefig(OUT + ".png", dpi=600)
fig.savefig(OUT + ".pdf")
plt.close()
# also print the top bands per motif for reference
for (ff, af), label in zip(FILES, LABELS):
    f, a = load(ff, af); amax = a.max()
    top = sorted([(round(fi,1), round(ai/amax*100,1)) for fi, ai in zip(f,a) if fi>50], key=lambda z:-z[1])[:6]
    print(f"{label}: {top}")
print(f"gas I2 = {fgas:.1f} cm-1 (exp 213)")
print(f"saved {OUT}.png/.pdf  (XLIM={XLIM})")
