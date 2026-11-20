#!/usr/bin/env python3
"""Final deliverables for the 442 Raman campaign: per-system mode tables and the
broadened-spectra figure vs the experimental bands (202-203 intact I2, ~167
chain-bridged I2). Read-only vs the campaign npy files; writes into
PBE_final_results_442/raman_442/.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv, os

C = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442'
OUT = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results_442/raman_442'
os.makedirs(OUT, exist_ok=True)
SYS = [('st_poly_lS5', 'I2@Y-MOFst (chains)', '0.03'),
       ('vc_poly_lS3', 'I2@Y-MOFvc (chains + I$_2$)', '0.02'),
       ('act_poly_lS2', 'I2@Y-MOFact (branched chains)', '0.04'),
       ('act_mol_rd', 'I2@Y-MOFact(mol.) (molecular I$_2$)', '0.04'),
       ('adduct', 'I$_2$ + carboxylate O (dissociation/re-attachment)', '0.04')]

# ---- mode table ----
rows = []
for tag, label, dg in SYS:
    f = np.load(f'{C}/{tag}/raman/freqs.npy')
    a = np.load(f'{C}/{tag}/raman/raman_activities.npy')
    order = np.argsort(f)[::-1]
    for k, ix in enumerate(order):
        if f[ix] < 5:
            continue
        rows.append((label, dg, k + 1, round(float(f[ix]), 2), round(float(a[ix]), 1)))
with open(f'{OUT}/raman_modes_442.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['system', 'smearing_Ry', 'mode#', 'freq_cm-1', 'rel_activity'])
    w.writerows(rows)
print(f'raman_modes_442.csv: {len(rows)} modes')

# ---- spectra figure ----
PANELS = [('st_poly_lS5', 'Y-MOF$_{st}$\n(chains)', '0.03'),
          ('vc_poly_lS3', 'Y-MOF$_{vc}$\n(chains + I$_2$)', '0.02'),
          ('act_poly_lS2', 'Y-MOF$_{act}$\n(branched chains)', '0.04'),
          ('act_mol_rd', 'Y-MOF$_{act}$(mol.)\n(molecular I$_2$)', '0.04'),
          ('adduct', 'carboxylate-O test\n(dissoc./re-att.)', '0.04')]
fig, axes = plt.subplots(len(PANELS), 1, figsize=(8.6, 12.2), sharex=True)
X = np.linspace(0, 240, 1200)
GAMMA = 6.0
for ax, (tag, label, dg) in zip(axes, PANELS):
    f = np.load(f'{C}/{tag}/raman/freqs.npy')
    a = np.load(f'{C}/{tag}/raman/raman_activities.npy')
    y = np.zeros_like(X)
    for fx, act_ in zip(f, a):
        if act_ < 0.5 or fx < 15:
            continue
        y += act_ * np.exp(-((X - fx) ** 2) / (2 * GAMMA ** 2))
    y *= 88.0 / y.max()          # normalize broadened curve: no frame overflow
    ax.fill_between(X, y, color='#2b6cb0', alpha=0.55, lw=0)
    ax.plot(X, y, color='#2b6cb0', lw=1.0)
    ax.text(-0.115, 0.5, label, transform=ax.transAxes, ha='center', va='center', fontsize=10)
    ax.set_ylim(0, 130)
    ax.set_yticks([])
    ax.axvspan(199, 214, color='#c05621', alpha=0.13, lw=0)
    ax.axvspan(160, 174, color='#38a169', alpha=0.13, lw=0)
axes[0].text(206.5, 118, 'I$_2$(mol.)', ha='center', fontsize=8, color='#c05621')
axes[0].text(167, 118, 'bridged I$_2$', ha='center', fontsize=8, color='#38a169')
axes[-1].set_xlabel('Raman shift (cm$^{-1}$)', fontsize=10)
axes[0].set_title('Computed I-sublattice Raman spectra of I$_2$@Y-MOF (FD Hessian, 4×4×2)\n'
                  'sticks broadened by 6 cm$^{-1}$; shaded guides: experimental bands',
                  fontsize=10.5)
fig.text(0.99, 0.005, 'smearing: st 0.03 / vc 0.02 / act, act(mol.) 0.04 Ry; '
         'relative activities, bond-polarizability model', ha='right', fontsize=7, color='0.35')
fig.subplots_adjust(left=0.19, right=0.985, top=0.94, bottom=0.05, hspace=0.18)
fig.savefig(f'{OUT}/raman_spectra_442.png', dpi=220)
fig.savefig(f'{OUT}/raman_spectra_442.pdf')
print('raman_spectra_442.png/.pdf written')
