import sys, os, numpy as np, ase.io
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
     "I":"I.pbe-n-kjpaw_psl.1.0.0.UPF","Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF",
     "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF","Cd":"Cd.pbe-n-kjpaw_psl.1.0.0.UPF","N":"N.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"I":126.904,"Gd":157.25,"Y":88.906,"Tb":158.925,"Cd":112.414,"N":14.007}
poscar,outdir,label=sys.argv[1],sys.argv[2],sys.argv[3]
a=ase.io.read(poscar); cell=np.array(a.cell); frac=a.get_scaled_positions(wrap=False); sy=a.get_chemical_symbols()
uniq=[]; [uniq.append(s) for s in sy if s not in uniq]
os.makedirs(outdir,exist_ok=True)
L=["&CONTROL","  calculation = 'scf'",f"  prefix = '{label}'","  outdir = './tmp'",
   "  pseudo_dir = '/workspace/qe-7.5/pseudo'","  tprnfor = .true.","/",
   "&SYSTEM","  ibrav = 0",f"  nat = {len(sy)}",f"  ntyp = {len(uniq)}","  ecutwfc = 50","  ecutrho = 400",
   "  occupations = 'fixed'","  nspin = 1",
   "  input_dft = 'hse'","  exxdiv_treatment = 'gygi-baldereschi'","  x_gamma_extrapolation = .true.","/",
   "&ELECTRONS","  conv_thr = 1.0d-7","  diagonalization = 'david'","  diago_david_ndim = 4",
   "  mixing_mode = 'local-TF'","  mixing_beta = 0.2","  mixing_ndim = 16","  electron_maxstep = 500","/",
   "ATOMIC_SPECIES"]
for s in uniq: L.append(f"{s:<3s} {MASS[s]:>10.3f}  {PSE[s]}")
L.append("CELL_PARAMETERS angstrom")
for r in cell: L.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}")
L.append("ATOMIC_POSITIONS crystal")
for s,f in zip(sy,frac): L.append(f"{s:<3s} {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}")
L+=["K_POINTS gamma",""]
open(f"{outdir}/scf.in","w").write("\n".join(L))
print(f"{label}: {len(sy)} atoms -> {outdir}/scf.in")
