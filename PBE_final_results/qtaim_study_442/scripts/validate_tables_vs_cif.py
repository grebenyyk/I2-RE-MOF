#!/usr/bin/env python3
"""Independent validation of the 442 BCP tables against the CIF geometries."""
import numpy as np, csv, re
T='/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/qtaim_study_442/tables'
CIFS={'Y-MOFst(poly)':'l442_st_poly_lS5_relaxed.cif','Y-MOFvc(poly)':'l442_vc_poly_lS3_relaxed.cif',
      'Y-MOFact(poly)':'l442_act_poly_lS2_relaxed.cif','Y-MOFact(mol)':'l442_act_mol_rd_relaxed.cif'}
CELLS={'Y-MOFst(poly)':np.array([[8.9348,0,0],[0,9.6608,0],[-3.7625602248,0,16.5958256665]]),
       'Y-MOFact(poly)':np.array([[8.9348,0,0],[0,9.6608,0],[-3.7625602248,0,16.5958256665]]),
       'Y-MOFact(mol)':np.array([[8.9348,0,0],[0,9.6608,0],[-3.7625602248,0,16.5958256665]]),
       'Y-MOFvc(poly)':np.array([[8.7975,0,0],[-1.4362,8.2029,0],[-4.1207,-0.9336,17.2148]])}
def load_cif(p):
    s=open(p).read(); lab=s.index('_atom_site_fract_z'); out={}
    for l in s[lab:].splitlines():
        q=l.split()
        if len(q)>=5 and q[1] in ('I','O'):
            out.setdefault(q[1],[]).append([float(q[2]),float(q[3]),float(q[4])])
    return out
def dist(cell,a,b):
    d=np.array(a)-np.array(b); d-=np.round(d)
    return np.linalg.norm(d@cell)
bad=0; n=0
for r in csv.DictReader(open(f'{T}/BCP_II_full_442.csv')):
    if r.get('System') not in CELLS: continue
    n+=1
    cell=CELLS[r['System']]
    cif=load_cif(f'/workspace/dimitry/calc/qe/i2_mof/PBE_final_results_442/landscape_3variants/cifs/{CIFS[r["System"]]}')
    m=re.match(r'I(\d+)-I(\d+)',r['Bond'])
    d=dist(cell,cif['I'][int(m.group(1))-1],cif['I'][int(m.group(2))-1])
    if abs(d-float(r['d_pair(A)']))>=0.01:
        bad+=1
        print(f"II MISMATCH {r['System']} {r['Bond']}: table {r['d_pair(A)']} vs CIF {d:.3f}")
print(f"BCP_II_full_442: {n} rows, {n-bad} match CIF exactly, {bad} mismatches")
bad=0; n=0
for r in csv.DictReader(open(f'{T}/BCP_IO_full_442.csv')):
    if r.get('System') not in CELLS: continue
    n+=1
    cell=CELLS[r['System']]
    cif=load_cif(f'/workspace/dimitry/calc/qe/i2_mof/PBE_final_results_442/landscape_3variants/cifs/{CIFS[r["System"]]}')
    d=float(r['d(A)'])
    if not any(abs(d-dist(cell,o,ii))<0.01 for o in cif['O'] for ii in cif['I']):
        bad+=1
        print(f"IO MISMATCH {r['System']} {r['Bond']}: d={d} not found in CIF O-I distances")
print(f"BCP_IO_full_442: {n} rows, {n-bad} found in CIF geometry, {bad} mismatches")
print("\n=== corrected vc I-I rows (d_pair / rho_BCP / BONDED / path_ok / src) ===")
for r in csv.DictReader(open(f'{T}/BCP_II_full_442.csv')):
    if r.get('System')=='Y-MOFvc(poly)':
        print(f"  {r['Bond']}: {r['d_pair(A)']} A  rho={r['rho_BCP(e/bohr3)']}  BONDED={r['BONDED(0.02)']}  path_ok={r['path_ok']}  ({r['src']})")
