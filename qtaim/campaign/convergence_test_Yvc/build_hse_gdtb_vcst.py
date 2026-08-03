#!/usr/bin/env python3
"""HSE06 for the Gd/Tb metallic variants (vc + st), same geometries as the PBE table.
Loaded (metallic) -> HSE 4x4x2 (expensive); empty (insulating) -> HSE Gamma (cheap).
vc loaded={m}_deep/scf_loaded.in, vc empty={m}_deep/vc_empty.out (relaxed, verified ref);
st loaded={st_dir}/uniform_chain/loaded.out, st empty={st_dir}/empty_mof.out.
gas-I2 HSE done. Mirrors build_hse_yvcst.py for Gd/Tb."""
import re,os,numpy as np
BASE="/workspace/dimitry/calc/qe/i2_mof"; CAMP=f"{BASE}/qtaim/campaign"; OUT=f"{CAMP}/convergence_test_Yvc"
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
     "Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF","Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Gd":157.25,"Tb":158.925,"I":126.904}
def cell_from_in(fn):
    L=open(fn).read().splitlines(); ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
    return np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)])
def geom_out_final(f):
    t=open(f).read(); hdrs=list(re.finditer(r"CELL_PARAMETERS",t)); cell=None
    if hdrs:
        rows=[]
        for l in t[hdrs[-1].start():].splitlines()[1:]:
            p=l.split()
            if len(p)>=3:
                try: rows.append([float(p[0]),float(p[1]),float(p[2])])
                except ValueError: pass
            if len(rows)==3: break
        if len(rows)==3: cell=np.array(rows)
    if cell is None: cell=cell_from_in(f.replace(".out",".in"))
    m=re.search(r"Begin final coordinates.*?ATOMIC_POSITIONS\s*\(crystal\)\s*\n(.*?)End final coordinates",t,re.DOTALL); pos=[]
    if m:
        for l in m.group(1).splitlines():
            s=l.split()
            if len(s)>=4 and s[0] in PSE:
                try: pos.append((s[0],float(s[1]),float(s[2]),float(s[3])))
                except: pass
    if not pos: return geom_inp(f.replace(".out",".in"))
    return cell,pos
def geom_inp(f):
    L=open(f).read().splitlines()
    ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
    ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][0]
    cell=np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)]); pos=[];k=ai+1
    while k<len(L):
        t=L[k].split()
        if not t or t[0] in ("K_POINTS","ATOMIC_SPECIES","CONSTRAINTS"): break
        if t[0] in PSE and len(t)>=4:
            try: pos.append((t[0],float(t[1]),float(t[2]),float(t[3])))
            except: pass
        k+=1
    return cell,pos
def emit(tag,cell,pos,loaded):
    uniq=[]; [uniq.append(e) for e in [p[0] for p in pos] if e not in uniq]
    if loaded: km="4 4 2 0 0 0"; kp="K_POINTS automatic"; dg="0.02"; mb="0.05"; nd=24
    else: km=""; kp="K_POINTS gamma"; dg="0.03"; mb="0.1"; nd=20
    L=[f"&CONTROL\n  calculation = 'scf'\n  prefix = '{tag}'\n  outdir = './tmp_{tag}'",
       "  pseudo_dir = '/workspace/qe-7.5/pseudo'\n  tprnfor = .true.\n/",
       f"&SYSTEM\n  ibrav = 0\n  nat = {len(pos)}\n  ntyp = {len(uniq)}\n  ecutwfc = 50\n  ecutrho = 400",
       f"  occupations = 'smearing'\n  smearing = 'mv'\n  degauss = {dg}\n  nspin = 1",
       "  input_dft = 'hse'\n  exxdiv_treatment = 'gygi-baldereschi'\n  x_gamma_extrapolation = .true.",
       "  exx_fraction = 0.25\n  screening_parameter = 0.106\n  ecutfock = 150",
       "  vdw_corr = 'grimme-d3'\n  dftd3_version = 4\n/",
       "&ELECTRONS\n  conv_thr = 1.0d-6\n  diagonalization = 'david'\n  diago_david_ndim = 2",
       f"  mixing_mode = 'local-TF'\n  mixing_beta = {mb}\n  mixing_ndim = {nd}\n  electron_maxstep = 200\n/",
       "ATOMIC_SPECIES"]
    for e in uniq: L.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell: L.append(f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for e,x,y,z in pos: L.append(f"  {e:<2s} {x:.10f} {y:.10f} {z:.10f}")
    L.append(f"{kp}"+(f"\n  {km}\n" if km else "\n"))
    open(f"{OUT}/{tag}.in","w").write("\n".join(L)); return tag
# metal: (tag, deep_dir, st_dir)
META=[("gd",f"{CAMP}/Gd_deep",f"{BASE}/08_Gd"),("tb",f"{CAMP}/Tb_deep",f"{BASE}/06_Tb")]
JOBS=[]
for tg,deep,std in META:
    JOBS+=[(f"hse_empty_{tg}vc",f"{deep}/vc_empty.out",False),
           (f"hse_empty_{tg}st",f"{std}/empty_mof.out",False),
           (f"hse_poly_{tg}vc", f"{deep}/scf_loaded.in",True),
           (f"hse_poly_{tg}st", f"{std}/uniform_chain/loaded.out",True)]
man=[]
for tag,src,loaded in JOBS:
    if not os.path.exists(src): print(f"[MISS] {tag}: {src}"); continue
    cell,pos = geom_out_final(src) if src.endswith(".out") else geom_inp(src)
    if cell is None or not pos: print(f"[FAIL parse] {tag}"); continue
    man.append(emit(tag,cell,pos,loaded))
    print(f"[ok] {tag:16s} nat={len(pos):3d} {'4x4x2' if loaded else 'Gamma':6s} src={src.split(BASE+'/')[1]}")
open(f"{OUT}/hse_gdtb_manifest.txt","w").write("\n".join(man)+"\n")
print(f"\nwrote {len(man)} inputs -> hse_gdtb_manifest.txt")
