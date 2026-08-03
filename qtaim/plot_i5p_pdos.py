#!/usr/bin/env python3
"""
Plot the summed I-5p PDOS for the I2@MOF top-3 (08_Gd, 06_Tb, 13_Y), shifted to E-E_F.
Reads QE projwfc.x output (pdos.pdos_atm#N(I)_wfc#2(p) files) and the Fermi energy.

Outputs:
  - i5p_pdos.png            : I 5p DOS vs E-E_F for the three systems
  - stdout                  : integrated I-5p occupation (occupied < E_F, per I)
                              and the empty I-5p (sigma*) weight 0..1.5 eV above E_F

Usage:  python plot_i5p_pdos.py
"""
import sys
import numpy as np
import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# system label, projwfc directory, Fermi energy (eV) from scf.out ("the Fermi energy is X ev")
SYSTEMS = [
    ("13_Y",  "qtaim/13Y_rho800_G", 2.1338),
    ("06_Tb", "qtaim/06Tb_qtaim",   2.0521),
    ("08_Gd", "qtaim/08Gd_qtaim",   1.9797),
]
COLORS = {"13_Y": "#1f77b4", "06_Tb": "#2ca02c", "08_Gd": "#d62728"}
OUTPNG = "qtaim/i5p_pdos.png"

# numpy>=2 renamed trapz -> trapezoid
trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


def load_i5p(directory):
    """Sum the I 5p (wfc#2(p)) ldos over all I atoms in the system.
    QE pdos file columns: # E(eV)  ldos(E)  pdos_m1  pdos_m2  pdos_m3
    ldos (col 2) = total p DOS for that atom."""
    files = sorted(glob.glob(directory + "/pdos.pdos_atm*(I)_wfc#2(p)"))
    tot, e, nI = None, None, 0
    for f in files:
        d = np.loadtxt(f, comments="#")
        e = d[:, 0]
        ldos = d[:, 1]
        tot = np.zeros_like(ldos) if tot is None else tot
        tot += ldos
        nI += 1
    return e, tot, nI


def main():
    fig, ax = plt.subplots(figsize=(9, 5.5))
    summary = []
    for name, dr, ef in SYSTEMS:
        e, tot, nI = load_i5p(dr)
        if nI == 0:
            print(f"WARNING: no I-5p pdos files found in {dr}", file=sys.stderr)
            continue
        esh = e - ef
        occ = trap(tot[esh <= 0], esh[esh <= 0])          # occupied I-5p (sigma + LP)
        mask = (esh > 0) & (esh < 1.5)
        above = trap(tot[mask], esh[mask])                 # empty I-5p (sigma*) near E_F
        summary.append((name, nI, occ, occ / nI, above))
        ax.plot(esh, tot, color=COLORS[name], lw=1.6, label=f"{name} ({nI} I atoms)")
    ax.axvline(0, color="k", ls="--", lw=1, label="E$_F$")
    ax.axvspan(-8, 0, color="0.85", zorder=0)
    ax.set_xlim(-8, 5)
    ax.set_ylim(0, None)
    ax.set_xlabel("E - E$_F$ (eV)")
    ax.set_ylabel("I 5p DOS (states/eV)")
    ax.set_title("I 5p PDOS (summed over all I) — top-3 pts trio")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPNG, dpi=130)
    print(f"saved {OUTPNG}")
    print(f"\n{'system':8s} {'I5p occ(<EF)':>13s} {'per I':>7s} {'I5p 0..1.5eV >EF (empty sigma*)':>32s}")
    for name, nI, occ, per, above in summary:
        print(f"{name:8s} {occ:13.3f} {per:7.3f} {above:32.3f}")


if __name__ == "__main__":
    main()
