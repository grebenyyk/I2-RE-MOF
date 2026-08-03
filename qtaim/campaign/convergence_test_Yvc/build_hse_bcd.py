#!/usr/bin/env python3
"""Phase B + C-light: HSE06/Gamma single-points for the two HSE tables.
  B  : Gd/Tb activated polyiodide (chain) + molecular            (4 calcs)
  C-l: empty framework x3 metals (insulating, cheap) + gas-I2    (4 calcs)
Y poly+mol already done (hse_poly_yact=-5767.38535, hse_mol_yact=-5767.22892).
All HSE06/Gamma, nspin=1 (f-in-core PPs), degauss 0.03 mv, 50/400, ecutfock 150.
Geometries = EXACT production relax-finals (sources from collect_final_table.py)."""
import re,os,numpy as np
BASE="/workspace/dimitry/calc/qe/i2_mof"; CAMP=f"{BASE}/qtaim/campaign"; OUT=f"{CAMP}/convergence_test_Yvc"
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
     "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF",
     "Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"Gd":157.25,"Tb":158.925,"I":126.904}
def cell_from_in(fn):
    L=open(fn).read().splitlines()
    ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
    return np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)])
def geom_out_final(f):
    """Final positions from 'Begin final coordinates'; cell from last CELL_PARAMETERS
    block in .out, else from companion .in (fixed-cell relax)."""
    t=open(f).read()
    cells=re.findall(r"CELL_PARAMETERS\s*\(angstrom\)\s*\n((?:\s*-?\d[\d.\s-]*\n){3})",t)
    cell=np.array([[float(x) for x in l.split()[:3]] for l in cells[-1].splitlines()]) if cells else cell_from_in(f.replace(".out",".in"))
    m=re.search(r"Begin final coordinates.*?ATOMIC_POSITIONS\s*\(crystal\)\s*\n(.*?)End final coordinates",t,re.DOTALL)
    pos=[]
    if m:
        for l in m.group(1).splitlines():
            s=l.split()
            if len(s)>=4 and s[0] in PSE:
                try: pos.append((s[0],float(s[1]),float(s[2]),float(s[3])))
                except: pass
    if not pos:  # relax unfinished -> use .in initial
        return geom_inp(f.replace(".out",".in"))
    return cell,pos
def geom_inp(f):
    L=open(f).read().splitlines()
    ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
    ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][0]
    cell=np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)])
    pos=[];k=ai+1
    while k<len(L):
        t=L[k].split()
        if not t or t[0] in ("K_POINTS","ATOMIC_SPECIES","CONSTRAINTS"): break
        if t[0] in PSE and len(t)>=4:
            try: pos.append((t[0],float(t[1]),float(t[2]),float(t[3])))
            except: pass
        k+=1
    return cell,pos
def emit(tag,cell,pos):
    uniq=[]
    for e in [p[0] for p in pos]:
        if e not in uniq: uniq.append(e)
    L=[f"&CONTROL\n  calculation = 'scf'\n  prefix = '{tag}'\n  outdir = './tmp_{tag}'",
       "  pseudo_dir = '/workspace/qe-7.5/pseudo'\n  tprnfor = .true.\n/",
       f"&SYSTEM\n  ibrav = 0\n  nat = {len(pos)}\n  ntyp = {len(uniq)}\n  ecutwfc = 50\n  ecutrho = 400",
       "  occupations = 'smearing'\n  smearing = 'mv'\n  degauss = 0.03\n  nspin = 1",
       "  input_dft = 'hse'\n  exxdiv_treatment = 'gygi-baldereschi'\n  x_gamma_extrapolation = .true.",
       "  exx_fraction = 0.25\n  screening_parameter = 0.106\n  ecutfock = 150",
       "  vdw_corr = 'grimme-d3'\n  dftd3_version = 4\n/",
       "&ELECTRONS\n  conv_thr = 1.0d-6\n  diagonalization = 'david'\n  diago_david_ndim = 2",
       "  mixing_mode = 'local-TF'\n  mixing_beta = 0.1\n  mixing_ndim = 20\n  electron_maxstep = 100\n/",
       "ATOMIC_SPECIES"]
    for e in uniq: L.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell: L.append(f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for e,x,y,z in pos: L.append(f"  {e:<2s} {x:.10f} {y:.10f} {z:.10f}")
    L.append("K_POINTS gamma\n")
    open(f"{OUT}/{tag}.in","w").write("\n".join(L))
    return tag
JOBS=[  # (tag, source)  -- sources from collect_final_table.py
    ("hse_poly_gdact", f"{BASE}/08_Gd/activated_chain/loaded.out"),   # Gd polyiodide (chain)
    ("hse_mol_gdact",  f"{BASE}/08_Gd/activated_i2/loaded.out"),     # Gd molecular
    ("hse_poly_tbact", f"{BASE}/06_Tb/activated_i2/loaded_ext.out"), # Tb polyiodide (chain)
    ("hse_mol_tbact",  f"{BASE}/06_Tb/activated_i2_mol/loaded.out"), # Tb molecular
    ("hse_empty_yact", f"{BASE}/13_Y/activate_empty/relax.out"),     # Y empty (insulating)
    ("hse_empty_gdact",f"{BASE}/08_Gd/activate_empty/relax.out"),    # Gd empty
    ("hse_empty_tbact",f"{BASE}/06_Tb/activate_empty/relax.out"),    # Tb empty
    ("hse_i2_vac",     f"{BASE}/i2_vac_ref/scf.in"),]                # gas I2 ref
man=[]
for tag,src in JOBS:
    if not os.path.exists(src): print(f"[MISS] {tag}: {src}"); continue
    cell,pos = geom_out_final(src) if src.endswith(".out") else geom_inp(src)
    if cell is None or not pos: print(f"[FAIL parse] {tag}: {src}"); continue
    man.append(emit(tag,cell,pos))
    ni=sum(1 for p in pos if p[0]=='I')
    print(f"[ok] {tag:16s} nat={len(pos):3d} I={ni}  src={src.split(BASE+'/')[1]}")
open(f"{OUT}/hse_bcd_manifest.txt","w").write("\n".join(man)+"\n")
print(f"\nwrote {len(man)} HSE06/Gamma inputs -> {OUT}/hse_bcd_manifest.txt")
