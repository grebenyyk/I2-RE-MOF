import sys,re,numpy as np
OUT="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF","Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"I":126.904}
rt,ht=sys.argv[1],sys.argv[2]
L=open(f"{OUT}/{rt}.in").read().splitlines()
ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
cell=np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)])
t=open(f"{OUT}/{rt}.out").read(); blk=t.split("ATOMIC_POSITIONS")[-1]; atoms=[]
for l in blk.splitlines():
    s=l.split()
    if len(s)>=4 and re.match(r'^[A-Z]',s[0]):
        try: atoms.append((s[0],float(s[1]),float(s[2]),float(s[3])))
        except: pass
uniq=[];[uniq.append(a[0]) for a in atoms if a[0] not in uniq]
LL=[f"&CONTROL\n  calculation='scf'\n  prefix='{ht}'\n  outdir='./tmp_{ht}'","  pseudo_dir='/workspace/qe-7.5/pseudo'\n  tprnfor=.true.\n/",
   f"&SYSTEM\n  ibrav=0\n  nat={len(atoms)}\n  ntyp={len(uniq)}\n  ecutwfc=50\n  ecutrho=400","  occupations='smearing'\n  smearing='mv'\n  degauss=0.02\n  nspin=1",
   "  input_dft='hse'\n  exxdiv_treatment='gygi-baldereschi'\n  x_gamma_extrapolation=.true.\n  exx_fraction=0.25\n  screening_parameter=0.106\n  ecutfock=100\n  nqx1=4\n  nqx2=4\n  nqx3=2\n  vdw_corr='grimme-d3'\n  dftd3_version=4\n/",
   "&ELECTRONS\n  conv_thr=1.0d-6\n  diagonalization='david'\n  diago_david_ndim=2\n  mixing_mode='local-TF'\n  mixing_beta=0.05\n  mixing_ndim=24\n  electron_maxstep=300\n/","ATOMIC_SPECIES"]
for e in uniq: LL.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
LL.append("CELL_PARAMETERS angstrom")
for r in cell: LL.append(f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
LL.append("ATOMIC_POSITIONS crystal")
for e,x,y,z in atoms: LL.append(f"  {e:<2s} {x:.10f} {y:.10f} {z:.10f}")
LL.append("K_POINTS automatic\n  4 4 2 0 0 0\n")
open(f"{OUT}/{ht}.in","w").write("\n".join(LL)); print(f"wrote {ht}.in ({len(atoms)} atoms)")
