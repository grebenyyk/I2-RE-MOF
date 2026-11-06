#!/usr/bin/env python3
"""SI table: anchoring vs local contacts — where the anchoring energy lives.

Reconciles the frozen-fragment decomposition (decomp_442) with the QTAIM
contact census (442 CP tables): |dE_anchor| vs the summed local energy
densities (EML -V/2) at the framework<->polyiodide BCPs, plus the D3
dispersion / electronic split of dEads. All energies kJ per I2 (4 I2/cell).
Sources are read, not hardcoded, so the table refreshes with the CSVs.
"""
import csv

DECOMP = '/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/decomp_442/decomposition_table_442.csv'
T442 = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/qtaim_study_442/tables'
ADS = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results_442/adsorption_9systems/adsorption_table_442.csv'
OUT = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/qtaim_study_442/tables'
N_I2 = 4  # I2 per cell (8 I atoms)

ORDER = [('st', 'st_poly_lS5', 'Y-MOFst(poly)', 'I2@Y-MOFst'),
         ('act', 'act_poly_lS2', 'Y-MOFact(poly)', 'I2@Y-MOFact'),
         ('vc', 'vc_poly_lS3', 'Y-MOFvc(poly)', 'I2@Y-MOFvc')]

def rows(f):
    return list(csv.DictReader(open(f)))

decomp = {(r['variant'], r['config']): r for r in rows(DECOMP)}
ads = {r['system']: r for r in rows(ADS)}

out = []
print(f"{'system':14s} {'dE_assm':>8s} {'dE_anch':>8s} {'dE_def':>7s} | "
      f"{'N_YI':>4s} {'E(YI)':>6s} {'N_IO':>4s} {'E(IO)':>6s} {'sumE_loc':>8s} | "
      f"{'nonlocal':>8s} {'share':>6s} | {'dEdisp':>7s} {'dEelec':>7s}")
for var, cfg, sys442, adsname in ORDER:
    d = decomp[(var, cfg)]
    assm, anch, df = float(d['dEassembly_kJ/I2']), float(d['dEanchor_kJ/I2']), float(d['dEdef_kJ/I2'])
    n_yi = sum(1 for r in rows(f'{T442}/BCP_MI_native_442.csv') if r['System'] == sys442)
    e_yi = sum(float(r['E=-V/2']) for r in rows(f'{T442}/BCP_MI_native_442.csv') if r['System'] == sys442)
    ic = rows(f'{T442}/BCP_IC_full_442.csv')
    n_ic = sum(1 for r in ic if r['System'] == sys442)
    e_ic = sum(float(r['E=-V/2']) for r in ic if r['System'] == sys442)
    n_io = sum(1 for r in rows(f'{T442}/BCP_IO_full_442.csv') if r['System'] == sys442)
    e_io = sum(float(r['E=-V/2']) for r in rows(f'{T442}/BCP_IO_full_442.csv') if r['System'] == sys442)
    n_int, e_int = n_yi + n_ic + n_io, e_yi + e_ic + e_io
    loc = e_int / N_I2                                   # per I2
    anch_nl = abs(anch) - loc
    share = anch_nl / abs(anch) * 100
    a = ads[adsname]
    disp, elec = float(a['dEdisp_I2']), float(a['delec_I2'])
    out.append(dict(variant=var, config=cfg, dEassm=assm, dEanch=anch, dEdef=df,
                    n_yi=n_yi, e_yi=e_yi, n_io=n_io, e_io=e_io, n_ic=n_ic, e_ic=e_ic,
                    n_int=n_int, loc=loc, nonlocal_anchor=anch_nl, share=share,
                    dEdisp=disp, dEelec=elec))
    print(f"{adsname:14s} {assm:8.2f} {anch:8.2f} {df:7.2f} | "
          f"{n_yi:4d} {e_yi:6.1f} {n_io:4d} {e_io:6.1f} {loc:8.2f} | "
          f"{anch_nl:8.2f} {share:5.1f}% | {disp:7.2f} {elec:7.2f}")

with open(f'{OUT}/anchoring_reconciliation_442.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['system', 'dE_assembly_kJ/I2', 'dE_anchor_kJ/I2', 'dE_def_kJ/I2',
                'N_YI', 'sumE_YI_kJ/cell', 'N_IO', 'sumE_IO_kJ/cell', 'N_CI', 'sumE_CI_kJ/cell',
                'sumE_interface_kJ/I2', 'nonlocal_anchor_kJ/I2', 'nonlocal_share_%',
                'dEdisp_kJ/I2', 'dEelec_kJ/I2'])
    for r in out:
        w.writerow([f"I2@Y-MOF{r['variant']}", f"{r['dEassm']:.2f}", f"{r['dEanch']:.2f}", f"{r['dEdef']:.2f}",
                    r['n_yi'], f"{r['e_yi']:.2f}", r['n_io'], f"{r['e_io']:.2f}",
                    r['n_ic'], f"{r['e_ic']:.2f}", f"{r['loc']:.2f}",
                    f"{r['nonlocal_anchor']:.2f}", f"{r['share']:.1f}",
                    f"{r['dEdisp']:.2f}", f"{r['dEelec']:.2f}"])
print(f"\nwrote {OUT}/anchoring_reconciliation_442.csv")

try:
    from docx import Document
    doc = Document()
    doc.add_heading('Table S. Anchoring of the polyiodide chains: fragment decomposition '
                    'vs local QTAIM contact energies (4x4x2, PBE-D3(BJ), kJ per I2)', level=2)
    p = doc.add_paragraph()
    p.add_run('dE_assembly / dE_anchor / dE_def from the frozen-fragment decomposition at the '
              'relaxed winner geometries. sumE columns: Espinosa-Molins-Lecomte local energies '
              '(-V/2, Abramov G/V from rho_BCP and lap_BCP) summed over all framework-polyiodide '
              'bond critical points in the cell (Y...I from the QTAIM BCP analysis; I...O; C...I), '
              'divided by the 4 I2 per cell - an upper bound for closed-shell contacts. '
              '"Non-local anchoring" = |dE_anchor| - sumE_interface: the share carried by terms '
              'with no topological signature in rho(r) (long-range Coulomb of the delocalized '
              'polarization/charge transfer + interface dispersion). dEdisp / dEelec decompose '
              'the total dEads (PBE-D3 difference and remainder) of the 9-system table.')
    t = doc.add_table(rows=1, cols=10)
    t.style = 'Table Grid'
    hdr = ['System', 'dE_assembly', 'dE_anchor', 'dE_def', 'Y...I: N, sum E', 'I...O: N, sum E',
           'sum E_interface', 'Non-local anchoring', 'share %', 'dEdisp / dEelec']
    for c, h in zip(t.rows[0].cells, hdr):
        c.text = h
    for r in out:
        c = t.add_row().cells
        c[0].text = f"I2@Y-MOF{r['variant']}"
        c[1].text = f"{r['dEassm']:.2f}"
        c[2].text = f"{r['dEanch']:.2f}"
        c[3].text = f"{r['dEdef']:.2f}"
        c[4].text = f"{r['n_yi']}, {r['e_yi']:.1f}"
        c[5].text = f"{r['n_io']}, {r['e_io']:.1f}"
        c[6].text = f"{r['loc']:.2f}"
        c[7].text = f"{r['nonlocal_anchor']:.2f}"
        c[8].text = f"{r['share']:.1f}"
        c[9].text = f"{r['dEdisp']:.2f} / {r['dEelec']:.2f}"
    doc.save(f'{OUT}/anchoring_reconciliation_442.docx')
    print(f"wrote {OUT}/anchoring_reconciliation_442.docx")
except ImportError:
    print('python-docx not available; CSV only')
