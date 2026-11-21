#!/usr/bin/env python3
"""Pair-decomposition analysis (user request 3a): decompose the delocalized
iodine modes into per-pair quantities and correlate with the QTAIM data.

For every bonded I-I pair (r <= 3.25 A) of every system:
  K_p  = internal-stretch force constant from the saved FD Hessian block
         (diatomic approximation: both atoms displaced along +/- the bond axis)
  nu_local = sqrt(K_p/mu)*RY_TO_CM, mu = m_I/2 (isolated-diatomic approx)
  nu_PED   = PED-weighted average of the observed bands containing this pair
  rho_BCP, E(-V/2) joined from the production BCP_II_full_442 table
Output: pairs_analysis_442.csv + correlation figure (nu_local vs rho_BCP,
nu_local vs -V/2, nu_PED vs nu_local) into PBE_final_results_442/raman_442/.
"""
import numpy as np
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv, os

C = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/raman_442'
TAB = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/qtaim_study_442/tables'
OUT = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results_442/raman_442'
RY_TO_CM = 3634.77
M_I = 126.904
MU = M_I / 2.0

SYS = [('st_poly_lS5', 'Y-MOFst(poly)'), ('vc_poly_lS3', 'Y-MOFvc(poly)'),
       ('act_poly_lS2', 'Y-MOFact(poly)'), ('act_mol_rd', 'Y-MOFact(mol)'),
       ('adduct', 'Y-MOFadduct(test)')]

# production QTAIM table: (system-name, bond-set, d, rho_BCP e/bohr3, E=-V/2 kJ)
qrows = {}
for r in list(csv.DictReader(open(f'{TAB}/BCP_II_full_442.csv'))):
    key = (r['System'], frozenset(re.findall(r'I(\d+)', r['Bond'])))
    qrows.setdefault(key, []).append((float(r['d_pair(A)']), float(r['rho_BCP(e/bohr3)']), float(r['E=-V/2'])))
import re

out_rows = []
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
colors = {'Y-MOFst(poly)': '#1f77b4', 'Y-MOFvc(poly)': '#d62728', 'Y-MOFact(poly)': '#2ca02c',
          'Y-MOFact(mol)': '#7f7f7f', 'Y-MOFadduct(test)': '#9467bd'}
for tag, sname in SYS:
    f = np.load(f'{C}/{tag}/raman/freqs.npy')
    a = np.load(f'{C}/{tag}/raman/raman_activities.npy')
    evecs = np.load(f'{C}/{tag}/raman/evecs.npy')
    cell = np.load(f'{C}/{tag}/raman/cell.npy')
    Ifrac = np.load(f'{C}/{tag}/raman/I_frac_ref.npy')
    Phi = np.load(f'{C}/{tag}/raman/Phi.npy')          # 24x24, Ry/bohr^2
    bonds = []
    for i in range(8):
        for j in range(i+1, 8):
            df = Ifrac[i]-Ifrac[j]; df -= np.round(df)
            r = np.linalg.norm(df@cell)
            if r <= 3.25:
                bonds.append((i, j, r, (df@cell)/r))
    for (i, j, r, u) in bonds:
        # internal-stretch force constant: pattern v (i: -u/2, j: +u/2)
        v = np.zeros(24)
        v[3*i:3*i+3] = -u/2
        v[3*j:3*j+3] = +u/2
        K = v @ Phi @ v                                # Ry/bohr^2
        nu_loc = np.sqrt(K/MU) * RY_TO_CM if K > 0 else -np.sqrt(-K/MU)*RY_TO_CM
        # PED over observed modes
        ped_num, ped_den = [], []
        contrib = []
        for k, ix in enumerate(np.argsort(f)[::-1]):
            if a[ix] < 1.0 or f[ix] < 15:
                continue
            e = evecs[:, ix]
            s = np.dot(e[3*j:3*j+3]-e[3*i:3*i+3], u)/np.sqrt(2)
            tot = sum(abs(np.dot(e[3*jj:3*jj+3]-e[3*ii:3*ii+3], uu)/np.sqrt(2))**2
                      for ii, jj, rr, uu in bonds)
            ped = s**2/tot if tot > 0 else 0.0
            ped_num.append(ped*f[ix]); ped_den.append(ped)
            if ped >= 0.10:
                contrib.append(f'{f[ix]:.0f}({ped*100:.0f}%)')
        nu_ped = sum(ped_num)/sum(ped_den) if sum(ped_den) > 0 else np.nan
        labs = f'I{i+1}-I{j+1}'
        q = qrows.get((sname, frozenset({str(i+1), str(j+1)})))
        if q:
            dd = sum(x[0] for x in q)/len(q)
            rho = sum(x[1] for x in q)/len(q)
            ev = sum(x[2] for x in q)/len(q)
        else:
            dd = rho = ev = float('nan')
        out_rows.append([sname, labs, f'{r:.3f}', f'{dd:.3f}', f'{rho:.4f}', f'{ev:.1f}',
                         f'{K:.4f}', f'{nu_loc:.1f}', f'{nu_ped:.1f}' if np.isfinite(nu_ped) else '',
                         '; '.join(contrib) if contrib else '-'])
        col = colors[sname]
        mk = 'o' if np.isfinite(rho) else 'x'
        if np.isfinite(rho):
            axes[0].scatter(rho, nu_loc, s=28, color=col, alpha=0.85, edgecolors='k', lw=.3, label=None)
            axes[1].scatter(ev, nu_loc, s=28, color=col, alpha=0.85, edgecolors='k', lw=.3)
        if np.isfinite(nu_ped):
            axes[2].scatter(nu_loc, nu_ped, s=28, color=col, alpha=0.85, edgecolors='k', lw=.3)
        axes[0].scatter([], [], color=col, label=sname.replace('Y-MOF', ''))  # legend proxy
for ax, xl, yl in zip(axes, ('rho_BCP (e/Bohr$^{-3}$)', 'E(I-I) = -V/2 (kJ/mol)', 'local $\\nu_p$ (cm$^{-1}$)'),
                      ('local $\\nu_p$ (cm$^{-1}$)', 'local $\\nu_p$ (cm$^{-1}$)', 'PED-avg band $\\nu$ (cm$^{-1}$)')):
    ax.set_xlabel(xl, fontsize=9); ax.set_ylabel(yl, fontsize=9)
    ax.tick_params(labelsize=8)
axes[0].legend(fontsize=7, frameon=False)
axes[2].plot(0, 0)
lims = [50, 240]
axes[2].plot(lims, lims, 'k--', lw=0.7)
axes[2].set_xlim(lims); axes[2].set_ylim(lims)
axes[2].text(60, 210, 'y=x', fontsize=8, color='0.4')
fig.suptitle('Per-pair decomposition: local diatomic frequencies vs QTAIM descriptors', fontsize=10.5)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(f'{OUT}/pairs_correlation_442.png', dpi=200)
fig.savefig(f'{OUT}/pairs_correlation_442.pdf')

with open(f'{OUT}/pairs_analysis_442.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['system', 'pair', 'd_FD_geom(A)', 'd_S11(A)', 'rho_BCP(e/Bohr3)', 'E=-V/2(kJ/mol)',
                'K_local(Ry/Bohr2)', 'nu_local(cm-1)', 'nu_PEDavg(cm-1)', 'PED>=10%: band(PED%)'])
    for r in out_rows:
        w.writerow(r)
print(f'pairs_analysis_442.csv: {len(out_rows)} pairs')
print(axes[0].get_legend_handles_labels()[1])

# act vs vc discriminator printout
for r in out_rows:
    if (r[0], r[1]) in (('Y-MOFvc(poly)', 'I1-I5'), ('Y-MOFact(poly)', 'I3-I4'),
                        ('Y-MOFact(poly)', 'I1-I2'), ('Y-MOFvc(poly)', 'I2-I3')):
        print(r[:2], 'nu_local =', r[7], 'nu_PED =', r[8])
