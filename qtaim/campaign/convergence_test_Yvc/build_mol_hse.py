#!/usr/bin/env python3
"""Once PBE molecular-seed relaxes (mseed_Y{v}_S{s}) finish, build HSE/4x4x2
single-point inputs (hmol_Y{v}_S{s}) from each relaxed molecular geometry.
Locked settings: 4x4x2, dg0.02, ef100, 50/400, nqx=4x4x2, HSE06+D3."""
import re,sys,numpy as np
OUT="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF","Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"I":126.904}
def geom_relax_final(f):
    """framework from .in (fixed-cell relax), final I+framework pos from last ATOMIC_POSITIONS of .out."""
    inL=open(f.replace(".out",".in")).read().splitlines()
    ci=[i for i,l in enumerate(inL) if l.startswith("CELL_PARAMETERS")][0]
    cell=np.array([[float(x) for x in inL[ci+1+k].split()[:3]] for k in range(3)])
    t=open(f).read(); blk=t.split("ATOMIC_POSITIONS")[-1]; atoms=[]
    for l in blk.splitlines():
        s=l.split()
        if len(s)>=4 and re.match(r'^[A-Z]',s[0]):
            try: atoms.append((s[0],float(s[1]),float(s[2]),float(s[3])))
            except: pass
    return cell,atoms
def emit(tag,cell,atoms):
    uniq=[]; [uniq.append(e) for e in [a[0] for a in atoms] if e not in uniq]
    L=[f"&CONTROL\n  calculation='scf'\n  prefix='{tag}'\n  outdir='./tmp_{tag}'",
       "  pseudo_dir='/workspace/qe-7.5/pseudo'\n  tprnfor=.true.\n/",
       f"&SYSTEM\n  ibrav=0\n  nat={len(atoms)}\n  ntyp={len(uniq)}\n  ecutwfc=50\n  ecutrho=400",
       "  occupations='smearing'\n  smearing='mv'\n  degauss=0.02\n  nspin=1",
       "  input_dft='hse'\n  exxdiv_treatment='gygi-baldereschi'\n  x_gamma_extrapolation=.true.",
       "  exx_fraction=0.25\n  screening_parameter=0.106\n  ecutfock=100\n  nqx1=4\n  nqx2=4\n  nqx3=2",
       "  vdw_corr='grimme-d3'\n  dftd3_version=4\n/",
       "&ELECTRONS\n  conv_thr=1.0d-6\n  diagonalization='david'\n  diago_david_ndim=2",
       "  mixing_mode='local-TF'\n  mixing_beta=0.05\n  mixing_ndim=24\n  electron_maxstep=300\n/",
       "ATOMIC_SPECIES"]
    for e in uniq: L.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell: L.append(f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for e,x,y,z in atoms: L.append(f"  {e:<2s} {x:.10f} {y:.10f} {z:.10f}")
    L.append("K_POINTS automatic\n  4 4 2 0 0 0\n")
    open(f"{OUT}/{tag}.in","w").write("\n".join(L))
variants=sys.argv[1:] if len(sys.argv)>1 else ["st","vc"]
man=[]
for v in variants:
    for s in range(1,2):
        src=f"{OUT}/mseed_Y{v}_S{s}.out"
        tag=f"hmol_Y{v}_S{s}"; cell,atoms=geom_relax_final(src)
        emit(tag,cell,atoms); man.append(tag)
        print(f"[ok] {tag}: {len(atoms)} atoms")
open(f"{OUT}/hmol_manifest.txt","w").write("\n".join(man)+"\n")
print(f"\nwrote {len(man)} HSE molecular inputs")
