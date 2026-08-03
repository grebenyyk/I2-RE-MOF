#!/usr/bin/env python3
"""
Per-system PDOS: I 5p (summed over I atoms) + metal d (Y 4d / Tb 5d / Gd 5d, summed
over metal atoms), vs E-E_F. Shows the I-5p <-> metal-d energy overlap below E_F
(the I(sigma) -> metal(d) sigma-donation mixing) and the empty sigma* above E_F.

Reads QE projwfc.x output. Channels: I wfc#2(p) = 5p; metal wfc#5(d) = d.
Usage: python plot_i5p_metald_pdos.py
"""
import numpy as np, glob
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

SYSTEMS = [("13_Y","Y","qtaim/13Y_rho800_G",2.1338),
           ("06_Tb","Tb","qtaim/06Tb_qtaim",2.0521),
           ("08_Gd","Gd","qtaim/08Gd_qtaim",1.9797)]
MCOLOR = {"Y":"#1f77b4","Tb":"#2ca02c","Gd":"#d62728"}
trap = np.trapezoid if hasattr(np,"trapezoid") else np.trapz

def sum_channel(dr, atom, wfc):
    """Sum ldos (col 2) over all atoms of species `atom` for the given wfc tag, e.g. '2(p)'."""
    files = sorted(glob.glob(f"{dr}/pdos.pdos_atm*({atom})_wfc#{wfc}"))
    tot, e, n = None, None, 0
    for f in files:
        d = np.loadtxt(f, comments="#"); e = d[:,0]; ldos = d[:,1]
        tot = np.zeros_like(ldos) if tot is None else tot; tot += ldos; n += 1
    return e, tot, n

fig, axes = plt.subplots(1, 3, figsize=(15,5), sharey=False)
print(f"{'system':8s} {'I5p occ':>8s} {'metald occ':>11s} {'metald 0..1.5eV>EF':>20s}")
for ax,(name,metal,dr,ef) in zip(axes, SYSTEMS):
    ei, i5p, nI = sum_channel(dr, "I", "2(p)")
    em, md, nM = sum_channel(dr, metal, "5(d)")
    # align on the energy grid (they share the same projwfc grid)
    e = ei; i5p = i5p; md = md
    esh = e - ef
    ax.fill_between(esh, i5p, color="0.7", step="mid", label=f"I 5p ({nI} I)")
    ax.plot(esh, md, color=MCOLOR[metal], lw=1.8, label=f"{metal} d ({nM} {metal})")
    ax.axvline(0, color="k", ls="--", lw=1); ax.axvspan(-8,0, color="0.9", zorder=0)
    ax.set_xlim(-8,5); ax.set_ylim(0,None)
    ax.set_xlabel("E - E$_F$ (eV)"); ax.set_title(name)
    ax.legend(loc="upper right", fontsize=8)
    iocc = trap(i5p[esh<=0], esh[esh<=0])
    mocc = trap(md[esh<=0], esh[esh<=0])
    mabove = trap(md[(esh>0)&(esh<1.5)], esh[(esh>0)&(esh<1.5)])
    print(f"{name:8s} {iocc:8.3f} {mocc:11.3f} {mabove:20.3f}")
axes[0].set_ylabel("DOS (states/eV)")
fig.suptitle("I 5p (gray) vs metal d (color) PDOS — overlap below E$_F$ = I(σ)→metal(d) mixing", y=1.02)
fig.tight_layout(); fig.savefig("qtaim/i5p_metald_pdos.png", dpi=130, bbox_inches="tight")
print("\nsaved qtaim/i5p_metald_pdos.png")
