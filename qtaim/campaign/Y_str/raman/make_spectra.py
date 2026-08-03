#!/usr/bin/env python3
"""Build the predicted Raman spectrum comparison plot (Y_str vs Y_act) and a per-system
mode table. Reads freqs.npy, evecs.npy, raman_activities.npy, cell.npy, I_frac_ref.npy
from each <sys>/raman/ directory.

Outputs (into Y_str/raman/):
  - spectrum_comparison.png / .pdf   (broadened + stem, 0-250 cm-1)
  - spectrum_comparison_log.png      (log-scale intensity, to show weak bands)
  - modes_table.txt                  (24 modes per system: freq, activity, bond character)
  - signature_bands.txt              (top Raman-active bands per system)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
SYSES = ["Y_str","Y_act"]
N_I = 8
SIGMA = 4.0   # cm-1 Lorentzian/Gaussian broadening half-width for spectrum

def load_sys(sysname):
    r = f"{BASE}/{sysname}/raman"
    freqs = np.load(f"{r}/freqs.npy")
    evecs = np.load(f"{r}/evecs.npy")
    cell = np.load(f"{r}/cell.npy")
    Ifrac = np.load(f"{r}/I_frac_ref.npy")
    try:
        act = np.load(f"{r}/raman_activities.npy")
    except FileNotFoundError:
        act = np.ones_like(freqs)  # uniform if no Raman model run
    # bonds
    bonds=[]
    for i in range(N_I):
        for j in range(i+1,N_I):
            df=Ifrac[i]-Ifrac[j]; df-=np.round(df)
            r=np.linalg.norm(df@cell)
            if r<4.0: bonds.append((i,j,r,(df@cell)/r))
    bonds.sort(key=lambda b:b[2])
    return freqs, evecs, act, bonds

def bond_character_for_mode(e, bonds, topK=3):
    """Return list of (bond_label, signed_stretch) for the top-K bonds in this mode."""
    stretches=[]
    for (i,j,r,u) in bonds:
        di = e[3*i:3*i+3]; dj = e[3*j:3*j+3]
        s = np.dot(dj-di, u)/np.sqrt(2)
        stretches.append((abs(s), i, j, r, s))
    stretches.sort(reverse=True)
    return stretches[:topK]

def broaden(freqs, acts, grid, sigma=SIGMA):
    """Lorentzian broadening."""
    # Lorentzian: L(x) = (1/pi) * (gamma/2) / ((x-x0)^2 + (gamma/2)^2)
    # use gamma = 2*sigma (FWHM) so half-width-at-half-max = sigma
    gamma = sigma
    spec = np.zeros_like(grid)
    for f,a in zip(freqs, acts):
        if f <= 0: continue  # skip imaginary / translational
        spec += a * (gamma/2.0)**2 / ((grid-f)**2 + (gamma/2.0)**2)
    return spec

def main():
    data = {}
    for s in SYSES:
        try:
            data[s] = load_sys(s)
            print(f"loaded {s}: {len(data[s][0])} modes, {len(data[s][3])} bonds")
        except FileNotFoundError as e:
            print(f"skip {s}: {e}")
    if len(data) < 2:
        print("need both systems; aborting plot"); sys.exit(1)

    # build mode table
    out_lines = []
    for s in SYSES:
        freqs, evecs, act, bonds = data[s]
        order = np.argsort(freqs)[::-1]
        out_lines.append(f"\n=== {s}: 24 I-sublattice modes (sorted high->low freq) ===")
        out_lines.append(f"{'mode':>4} {'freq':>9} {'activity':>10}  top bonds (I-I dist, signed stretch)")
        for k, ix in enumerate(order):
            e = evecs[:, ix]
            top = bond_character_for_mode(e, bonds)
            desc = "  ".join([f"I{i+1}-I{j+1}({r:.2f}A){'+' if s>=0 else '-'}{abs(s):.2f}"
                              for _,i,j,r,s in top])
            flag = " [IMAG]" if freqs[ix] < 0 else (" [near0]" if abs(freqs[ix])<10 else "")
            out_lines.append(f"{k+1:>4} {freqs[ix]:>9.2f} {act[ix]:>10.2f}  {desc}{flag}")
        # bond list
        out_lines.append(f"  bonds <4A: " + ", ".join([f"I{i+1}-I{j+1}={r:.2f}" for i,j,r,_ in bonds]))
    table = "\n".join(out_lines)
    print(table)
    with open(f"{BASE}/Y_str/raman/modes_table.txt","w") as fh: fh.write(table)

    # signature bands: top-5 Raman active modes per system
    sig_lines=[]
    for s in SYSES:
        freqs, evecs, act, bonds = data[s]
        sig_lines.append(f"\n=== {s}: top Raman-active modes ===")
        # only real modes
        mask = freqs > 5
        ix_real = np.where(mask)[0]
        order = ix_real[np.argsort(act[ix_real])[::-1]]
        for k, ix in enumerate(order[:6]):
            e = evecs[:, ix]
            top = bond_character_for_mode(e, bonds)
            desc = "  ".join([f"I{i+1}-I{j+1}({r:.2f}A){'+' if s2>=0 else '-'}{abs(s2):.2f}"
                              for _,i,j,r,s2 in top])
            sig_lines.append(f"  {freqs[ix]:>7.2f} cm-1  act={act[ix]:>6.2f}  {desc}")
    sig = "\n".join(sig_lines)
    print(sig)
    with open(f"{BASE}/Y_str/raman/signature_bands.txt","w") as fh: fh.write(sig)

    # spectrum plot
    grid = np.linspace(0, 250, 2500)
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    colors = {"Y_str":"#1f77b4", "Y_act":"#d62728"}
    # top: broadened spectra
    ax = axes[0]
    for s in SYSES:
        freqs, evecs, act, bonds = data[s]
        spec = broaden(freqs, act, grid)
        if spec.max()>0: spec = spec/spec.max()  # normalize each to 1
        ax.plot(grid, spec, color=colors[s], lw=2, label=s)
        # stem the actual mode positions
        real = freqs > 5
        ax.stem(freqs[real], act[real]/act[real].max() if act[real].max()>0 else act[real],
                linefmt=colors[s]+':', markerfmt=colors[s]+'o', basefmt=' ',
                label=f"{s} modes")
    ax.set_ylabel("Raman intensity (rel.)")
    ax.set_title("Predicted Raman spectrum of I sublattice in I2@Y-MOF (framework frozen)\n"
                 "PBE/DNP, finite-difference phonons, bond-polarizability model")
    ax.legend(loc="upper right", fontsize=9)
    ax.set_xlim(0, 250); ax.grid(alpha=0.3)
    # annotate reference bands
    for x,lab in [(110,"I3- sym"),(145,"I3- asym"),(165,"I5-"),(213,"I2 gas")]:
        ax.axvline(x, color='gray', ls='--', alpha=0.4)
        ax.text(x, 0.95, lab, rotation=90, fontsize=7, color='gray', va='top', ha='right')
    # bottom: stem only (cleaner)
    ax = axes[1]
    for s in SYSES:
        freqs, evecs, act, bonds = data[s]
        real = freqs > 5
        markerline, stemlines, baseline = ax.stem(freqs[real],
                act[real], linefmt=colors[s], markerfmt=colors[s]+'o', basefmt=' ',
                label=s)
        plt.setp(stemlines, 'linewidth', 1.5)
    ax.set_xlabel("Raman shift (cm-1)")
    ax.set_ylabel("Raman activity (arb.)")
    ax.legend(loc="upper right")
    ax.set_xlim(0, 250); ax.grid(alpha=0.3)
    for x,lab in [(110,"I3- sym"),(145,"I3- asym"),(165,"I5-"),(213,"I2 gas")]:
        ax.axvline(x, color='gray', ls='--', alpha=0.4)
    plt.tight_layout()
    for ext in ("png","pdf"):
        fig.savefig(f"{BASE}/Y_str/raman/spectrum_comparison.{ext}", dpi=150)
    # also a log version
    fig2, ax = plt.subplots(1,1, figsize=(10,5))
    for s in SYSES:
        freqs, evecs, act, bonds = data[s]
        spec = broaden(freqs, act, grid)
        if spec.max()>0: spec = spec/spec.max()
        ax.semilogy(grid, spec+1e-3, color=colors[s], lw=2, label=s)
    ax.set_xlim(0,250); ax.set_ylim(1e-3, 2)
    ax.set_xlabel("Raman shift (cm-1)"); ax.set_ylabel("log Raman intensity (rel.)")
    ax.legend(); ax.grid(alpha=0.3, which='both')
    ax.set_title("Predicted Raman (log scale)")
    plt.tight_layout()
    fig2.savefig(f"{BASE}/Y_str/raman/spectrum_comparison_log.png", dpi=150)
    print("saved spectrum_comparison.{png,pdf} and _log.png")

if __name__ == "__main__":
    main()
