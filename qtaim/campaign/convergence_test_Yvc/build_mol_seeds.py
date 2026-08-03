#!/usr/bin/env python3
"""Molecular landscape seeds: PBE/4x4x2 relax of PERTURBED polyiodide@variant starts.
PBE dissociates polyiodide -> molecular; 5 perturbations sample distinct molecular
minima in the same pore. Framework frozen (if_pos 0 0 0), I free (1 1 1).
usage: build_mol_seeds.py <variant>   (st|act|vc)  -> mseed_Y{v}_S{1..5}.in"""
import sys,re,numpy as np
BASE="/workspace/dimitry/calc/qe/i2_mof"; OUT=f"{BASE}/qtaim/campaign/convergence_test_Yvc"
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF","Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"I":126.904}
v=sys.argv[1] if len(sys.argv)>1 else "st"
# read polyiodide@variant geometry (framework + 8 I chain)
L=open(f"{OUT}/rd_poly_Y{v}.in").read().splitlines()
ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
cell=np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)])
ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][0]; atoms=[];k=ai+1
while k<len(L):
    s=L[k].split()
    if not s or s[0] in("K_POINTS","ATOMIC_SPECIES"): break
    if len(s)>=4 and s[0] in PSE:
        try: atoms.append([s[0],float(s[1]),float(s[2]),float(s[3])])
        except: pass
    k+=1
nfw=len(atoms)-8  # framework atoms frozen, last 8 I free
def emit(tag,atoms,kick_seed):
    rng=np.random.default_rng(kick_seed)
    uniq=[]; [uniq.append(a[0]) for a in atoms if a[0] not in uniq]
    # perturb the 8 I (fractional) by ~+-0.02 (fractional ~ 0.2-0.3 A)
    a2=[a[:] for a in atoms]
    for j in range(nfw,len(atoms)):
        a2[j][1]=(a2[j][1]+rng.uniform(-0.02,0.02))%1.0
        a2[j][2]=(a2[j][2]+rng.uniform(-0.02,0.02))%1.0
        a2[j][3]=(a2[j][3]+rng.uniform(-0.02,0.02))%1.0
    L=[f"&CONTROL\n  calculation='relax'\n  prefix='{tag}'\n  outdir='./tmp_{tag}'",
       "  pseudo_dir='/workspace/qe-7.5/pseudo'\n  tprnfor=.true.\n  etot_conv_thr=1.0d-5\n  forc_conv_thr=1.0d-3\n/",
       f"&SYSTEM\n  ibrav=0\n  nat={len(atoms)}\n  ntyp={len(uniq)}\n  ecutwfc=50\n  ecutrho=400",
       "  occupations='smearing'\n  smearing='mv'\n  degauss=0.02\n  nspin=1\n  vdw_corr='grimme-d3'\n  dftd3_version=4\n/",
       "&ELECTRONS\n  conv_thr=1.0d-7\n  diagonalization='david'\n  diago_david_ndim=2",
       "  mixing_mode='local-TF'\n  mixing_beta=0.2\n  mixing_ndim=24\n  electron_maxstep=300\n/",
       "&IONS\n  ion_dynamics='bfgs'\n  ion_maxstep=100\n/","ATOMIC_SPECIES"]
    for e in uniq: L.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell: L.append(f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for j,a in enumerate(a2):
        flags="0 0 0" if j<nfw else "1 1 1"
        L.append(f"  {a[0]:<2s} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}  {flags}")
    L.append("K_POINTS automatic\n  4 4 2 0 0 0\n")
    open(f"{OUT}/{tag}.in","w").write("\n".join(L))
for s in range(1,6):
    emit(f"mseed_Y{v}_S{s}",atoms,kick_seed=100+s)
print(f"wrote 5 perturbed polyiodide@Y{v} molecular-seed relax inputs (mseed_Y{v}_S1..S5)")
print(f"  framework frozen (first {nfw}), I free (last 8), 4x4x2 PBE+D3")
