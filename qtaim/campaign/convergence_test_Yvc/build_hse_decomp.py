#!/usr/bin/env python3
"""HSE chain-form/anchoring/deformation decomposition (corrected Born-Haber cycle).
Needs E*_I (iodine sublattice at loaded geom, no framework) for Yst/Yact/Yvc,
and E*_host (framework at loaded geom, no I) for Yvc only (st/act: Host*=Host0=rd_empty).
Geometries extracted from rd_poly_*.in (the loaded polyiodide@variant).
E*_I: 4x4x2 (may be metallic chain); E*_host: Gamma (insulator, matches rd_empty)."""
import re,os,numpy as np
BASE="/workspace/dimitry/calc/qe/i2_mof"; OUT=f"{BASE}/qtaim/campaign/convergence_test_Yvc"
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
     "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF","Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"Gd":157.25,"Tb":158.925,"I":126.904}
def rdgeom(tag,only=None,exclude=None):
    L=open(f"{OUT}/{tag}.in").read().splitlines()
    ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
    cell=np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)])
    ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][0]; atoms=[];k=ai+1
    while k<len(L):
        s=L[k].split()
        if not s or s[0] in("K_POINTS","ATOMIC_SPECIES"): break
        if len(s)>=4 and s[0] in PSE:
            e=s[0]
            if (only and e not in only) or (exclude and e in exclude): k+=1; continue
            try: atoms.append((e,float(s[1]),float(s[2]),float(s[3])))
            except: pass
        k+=1
    return cell,atoms
def emit(tag,cell,atoms,kmesh,nqx,note=""):
    uniq=[]; [uniq.append(e) for e in [a[0] for a in atoms] if e not in uniq]
    L=[f"&CONTROL\n  calculation='scf'\n  prefix='{tag}'\n  outdir='./tmp_{tag}'",
       "  pseudo_dir='/workspace/qe-7.5/pseudo'\n  tprnfor=.true.\n/",
       f"&SYSTEM\n  ibrav=0\n  nat={len(atoms)}\n  ntyp={len(uniq)}\n  ecutwfc=50\n  ecutrho=400",
       "  occupations='smearing'\n  smearing='mv'\n  degauss=0.02\n  nspin=1",
       "  input_dft='hse'\n  exxdiv_treatment='gygi-baldereschi'\n  x_gamma_extrapolation=.true.",
       "  exx_fraction=0.25\n  screening_parameter=0.106\n  ecutfock=100",
       f"  nqx1={nqx[0]}\n  nqx2={nqx[1]}\n  nqx3={nqx[2]}",
       "  vdw_corr='grimme-d3'\n  dftd3_version=4\n/",
       "&ELECTRONS\n  conv_thr=1.0d-6\n  diagonalization='david'\n  diago_david_ndim=2",
       "  mixing_mode='local-TF'\n  mixing_beta=0.05\n  mixing_ndim=24\n  electron_maxstep=300\n/",
       "ATOMIC_SPECIES"]
    for e in uniq: L.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell: L.append(f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for e,x,y,z in atoms: L.append(f"  {e:<2s} {x:.10f} {y:.10f} {z:.10f}")
    L.append(f"K_POINTS automatic\n  {kmesh}\n")
    open(f"{OUT}/{tag}.in","w").write("\n".join(L))
    return f"{tag}: {len(atoms)} atoms {note}"
man=[]
# E*_I (iodine sublattice, loaded geom, no framework) for Yst/Yact/Yvc -- 4x4x2
for v in ["st","act","vc"]:
    cell,atoms=rdgeom(f"rd_poly_Y{v}",only={"I"})
    man.append(emit(f"dc_isubl_Y{v}",cell,atoms,"4 4 2 0 0 0",(4,4,2),"I-only (chain)"))
# E*_host (framework at loaded geom, no I) for Yvc only -- Gamma
cell,atoms=rdgeom("rd_poly_Yvc",exclude={"I"})
man.append(emit("dc_hostld_Yvc",cell,atoms,"1 1 1 0 0 0",(1,1,1),"framework@loaded, vc"))
for m in man: print("wrote",m)
open(f"{OUT}/decomp_manifest.txt","w").write("\n".join(["dc_isubl_Yst","dc_isubl_Yact","dc_isubl_Yvc","dc_hostld_Yvc"])+"\n")
