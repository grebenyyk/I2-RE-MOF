#!/usr/bin/env python3
"""Table-2-style topological summary of the 442 CP tables: per system, the
I-I and I...O BCPs as (d range / number of unique BCPs) and (rho range
(mean)), rho in e x A^-3 — the same layout as the user's Gamma-era Table 2.

Definitions: 'unique' = distinct atom-pair bonds in the cell (periodic images
not multiplied); I-I set = BCPs with rho_BCP >= 0.02 a.u. (the tables' BONDED
criterion, i.e. genuine polyiodide bonding incl. the vc borderline 3.40/3.43
bridges); I...O set = every I...O BCP located (no rho cutoff).
"""
import csv
BOHR3_TO_A3 = 1 / 0.5291772109**3          # e/bohr^3 -> e/A^3 (= 6.74833)
T = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/qtaim_study_442/tables'
OUT = '/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/qtaim_study_442/tables'

# II table uses label pairs; collapse rows that repeat the same bond (dedup by
# the pair, keeping both whenever values differ <2% is unnecessary — the dups
# in the table are the same CP reported by two footings; take the first).
def uniq_rows(fname, key, bonded_only):
    seen, out = set(), []
    for r in csv.DictReader(open(f'{T}/{fname}')):
        if bonded_only and r.get('BONDED(0.02)') == '0':
            continue
        k = (r['System'], r[key])
        if k in seen:
            continue
        seen.add(k)
        out.append(r)
    return out

def summ(rows):
    ds = [float(r.get('d_pair(A)', r.get('d(A)'))) for r in rows]
    rs = [float(r['rho(e/A3)']) for r in rows]
    if not rows:
        return '', ''
    dr = (f'{min(ds):.2f}' if round(min(ds), 2) == round(max(ds), 2)
          else f'{min(ds):.2f}-{max(ds):.2f}')
    return (f'{dr} / {len(rows)}',
            f'{min(rs):.3f}-{max(rs):.3f} ({sum(rs)/len(rs):.3f})')

SYSORDER = [('Y-MOFst(poly)',  'Solvent-templated'),
            ('Y-MOFact(poly)', 'Activated'),
            ('Y-MOFact(mol)',  'Activated (molecular I$_2$)'),
            ('Y-MOFvc(poly)',  'Variable-cell')]

ii = uniq_rows('BCP_II_full_442.csv', 'Bond', bonded_only=True)
io = uniq_rows('BCP_IO_full_442.csv', 'Bond', bonded_only=False)

rows_out = []
print(f"{'system':28s} {'I-I d/N':22s} {'I-I rho':26s} {'I-O d/N':22s} {'I-O rho':26s}")
for sysname, group in SYSORDER:
    a = summ([r for r in ii if r['System'] == sysname])
    b = summ([r for r in io if r['System'] == sysname])
    rows_out.append((group, sysname, *a, *b))
    print(f'{sysname:28s} {a[0]:22s} {a[1]:26s} {b[0]:22s} {b[1]:26s}')

with open(f'{OUT}/topo_summary_442.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Group', 'System',
                'I-I d(A) range / N unique', 'I-I rho(e/A3) range (mean)',
                'I-O d(A) range / N unique', 'I-O rho(e/A3) range (mean)'])
    w.writerows(rows_out)
print(f"\nwrote {OUT}/topo_summary_442.csv")

# docx mirroring the user's Table 2 layout (System | I-I d/N | I-I rho | I-O d/N | I-O rho)
try:
    from docx import Document
    doc = Document()
    doc.add_heading('Table. Topological properties at the I-I and I...O bond critical points '
                    'in I2@Y-MOF (4x4x2 k-mesh, PBE+D3)', level=2)
    t = doc.add_table(rows=2, cols=5)
    t.style = 'Table Grid'
    top = t.rows[0].cells
    top[0].text = 'System'
    top[1].merge(top[2]).text = 'I-I'
    top[3].merge(top[4]).text = 'I...O'
    sub = t.rows[1].cells
    sub[0].text = ''
    sub[1].text = 'd, A; range / N unique'
    sub[2].text = 'rho_BCP, e- x A-3; range (mean)'
    sub[3].text = 'd, A; range / N unique'
    sub[4].text = 'rho_BCP, e- x A-3; range (mean)'
    last_group = None
    for group, sysname, di, ri, do, ro in rows_out:
        if group != last_group:
            c = t.add_row().cells
            c[0].text = group
            last_group = group
        c = t.add_row().cells
        c[0].text = sysname.replace('(poly)', '').replace('(mol)', '(mol.)')
        c[1].text = di
        c[2].text = ri
        c[3].text = do
        c[4].text = ro
    doc.save(f'{OUT}/topo_summary_442.docx')
    print(f'wrote {OUT}/topo_summary_442.docx')
except ImportError:
    print('python-docx not available; CSV only')
