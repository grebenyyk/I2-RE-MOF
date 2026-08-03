#!/usr/bin/env python3
"""Extract the 7 image geometries from the best NEB path (gd_neb_run2.path, iter~60,
force ~0.7) and write standalone single-point SCF inputs for each. Single-points have
no force-feedback loop, so they don't suffer the Broyden/NaN cascade that killed NEB.
conv_thr=1e-4 (partial-chain images floor ~3.8e-5 < 1e-4 => converge cleanly)."""
import re, numpy as np, os
BASE = "/workspace/dimitry/calc/qe/i2_mof"
GD = f"{BASE}/gd_neb"
BOHR = 0.5291772107  # bohr -> angstrom
cell = np.array([[9.0752,0,0],[0,9.6427,0],[-3.8981733929,0,16.7845745817]])
inv = np.linalg.inv(cell)
# species order from pw_1.in (C x32, Gd x4, H x12, O x24, I x8)
def species_order():
    t=open(f"{GD}/pw_1.in.bak").read().split("ATOMIC_POSITIONS",1)[1]
    syms=[]
    for l in t.splitlines():
        m=re.match(r"^\s*([A-Z][a-z]?)\s",l)
        if m and m.group(1) in ("C","Gd","H","O","I"): syms.append(m.group(1))
    return syms
syms=species_order()
assert len(syms)==80, len(syms)
# parse .path: Image blocks, each = energy line + 80 position lines (x y z fxs, bohr)
t=open(f"{GD}/gd_neb_run2.path").read()
images={}
for m in re.finditer(r"Image:\s*(\d+)\n\s*([-0-9.]+)\n((?:\s*[-0-9.]+\s+[-0-9.]+\s+[-0-9.]+.*\n){80})", t):
    n=int(m.group(1)); E=float(m.group(2))
    pos=[]
    for l in m.group(3).splitlines():
        x,y,z=float(l.split()[0]),float(l.split()[1]),float(l.split()[2])
        pos.append([x*BOHR,y*BOHR,z*BOHR])  # bohr->angstrom
    images[n]=(E,np.array(pos))
assert sorted(images)==[1,2,3,4,5,6,7], images.keys()
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF",
     "O":"O.pbe-kjpaw.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"Gd":157.25,"H":1.008,"O":15.999,"I":126.904}
print("image : E_neb(Ha) -> E(Ry) -> E(eV)")
for n in range(1,8):
    E,pos=images[n]
    frac=pos@inv
    pre=[s for s in sorted(set(syms),key=list(syms).index)]
    cnt=[syms.count(s) for s in pre]
    L=[f"&CONTROL",f"  calculation='scf'",f"  prefix='sp_img{n}'",f"  outdir='{GD}/tmp_sp/sp_img{n}'",
       f"  pseudo_dir='/workspace/qe-7.5/pseudo'","  tprnfor=.false.","/",
       "&SYSTEM","  ibrav=0","  nat=80",f"  ntyp={len(pre)}","  ecutwfc=50","  ecutrho=400",
       "  occupations='smearing'","  smearing='mv'","  degauss=0.01","  nspin=1",
       "  vdw_corr='grimme-d3'","  dftd3_version=4","/",
       "&ELECTRONS","  conv_thr=1.0d-4","  diagonalization='david'","  diago_david_ndim=4",
       "  mixing_mode='local-TF'","  mixing_beta=0.1","  mixing_ndim=16","  electron_maxstep=200","/",
       "ATOMIC_SPECIES"]
    for s in pre: L.append(f"{s:<3s} {MASS[s]:>10.3f}  {PSE[s]}")
    L+=["CELL_PARAMETERS angstrom"]
    for r in cell: L.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}")
    L+=["ATOMIC_POSITIONS crystal"]
    for s,f in zip(syms,frac): L.append(f"{s:<3s} {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}")
    L+=["K_POINTS gamma",""]
    open(f"{GD}/scf_img{n}.in","w").write("\n".join(L))
    print(f"  {n} : {E:.6f} -> {E*2:.4f} Ry -> {E*2*13.605693:.3f} eV   (wrote scf_img{n}.in)")
os.makedirs(f"{GD}/tmp_sp",exist_ok=True)
print("done")
