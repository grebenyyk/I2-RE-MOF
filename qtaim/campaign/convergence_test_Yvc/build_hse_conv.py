#!/usr/bin/env python3
"""CIF -> QE HSE06 input converter + Y-vc convergence sweep generator.
Reads authoritative CIFs (final_activated_results), emits HSE inputs with
explicit k-mesh AND nqx1-3 (EXX q-mesh), ecutfock, ecutwfc/rho. Sweeps on Y-vc
loaded (most metallic) to lock uniform settings for the full redo."""
import re,math,numpy as np,os
BASE="/workspace/dimitry/calc/qe/i2_mof"; CIF=f"{BASE}/final_activated_results"; OUT=f"{BASE}/qtaim/campaign/convergence_test_Yvc"
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
     "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF",
     "Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"Gd":157.25,"Tb":158.925,"I":126.904}
def cif_geom(f):
    t=open(f).read()
    def g(p): return float(re.search(p,t).group(1))
    a,b,c=[g(rf"_cell_length_{n}\s+([\d.\-eE]+)") for n in("a","b","c")]
    al,be,ga=[math.radians(g(rf"_cell_angle_{n}\s+([\d.\-eE]+)")) for n in("alpha","beta","gamma")]
    cell=np.zeros((3,3)); cell[0]=[a,0,0]; cell[1]=[b*math.cos(ga),b*math.sin(ga),0]
    cg=(math.cos(al)-math.cos(be)*math.cos(ga))/math.sin(ga)
    cell[2]=[c*math.cos(be),c*cg,c*math.sqrt(max(0,1-math.cos(be)**2-cg**2))]
    # atoms: lines with label type fx fy fz
    atoms=[]
    for l in t.splitlines():
        s=l.split()
        if len(s)>=5 and s[1] in PSE:
            try: atoms.append((s[1],float(s[2]),float(s[3]),float(s[4])))
            except: pass
    return cell,atoms
def emit(tag,cell,atoms,kmesh,nqx,ecutfock,ew,er,dg):
    uniq=[]; [uniq.append(e) for e in [a[0] for a in atoms] if e not in uniq]
    L=[f"&CONTROL\n  calculation = 'scf'\n  prefix = '{tag}'\n  outdir = './tmp_{tag}'",
       "  pseudo_dir = '/workspace/qe-7.5/pseudo'\n  tprnfor = .true.\n/",
       f"&SYSTEM\n  ibrav = 0\n  nat = {len(atoms)}\n  ntyp = {len(uniq)}\n  ecutwfc = {ew}\n  ecutrho = {er}",
       f"  occupations = 'smearing'\n  smearing = 'mv'\n  degauss = {dg}\n  nspin = 1",
       "  input_dft = 'hse'\n  exxdiv_treatment = 'gygi-baldereschi'\n  x_gamma_extrapolation = .true.",
       f"  exx_fraction = 0.25\n  screening_parameter = 0.106\n  ecutfock = {ecutfock}",
       f"  nqx1 = {nqx[0]}\n  nqx2 = {nqx[1]}\n  nqx3 = {nqx[2]}",
       "  vdw_corr = 'grimme-d3'\n  dftd3_version = 4\n/",
       "&ELECTRONS\n  conv_thr = 1.0d-6\n  diagonalization = 'david'\n  diago_david_ndim = 2",
       "  mixing_mode = 'local-TF'\n  mixing_beta = 0.05\n  mixing_ndim = 24\n  electron_maxstep = 300\n/",
       "ATOMIC_SPECIES"]
    for e in uniq: L.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell: L.append(f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for e,x,y,z in atoms: L.append(f"  {e:<2s} {x:.10f} {y:.10f} {z:.10f}")
    L.append(f"K_POINTS automatic\n  {kmesh}\n")
    open(f"{OUT}/{tag}.in","w").write("\n".join(L))
    return tag
# Y-vc loaded (most metallic) sweep
cell,atoms=cif_geom(f"{CIF}/vc-relaxed/Y_vc_loaded.cif")
print(f"Y_vc_loaded: {len(atoms)} atoms, {sum(1 for a in atoms if a[0]=='I')} I")
# (tag, k, nqx, ecutfock, ew, er)
SWEEP=[
 ("conv_yvc_kG",  "1 1 1 0 0 0",(1,1,1),150,50,400),
 ("conv_yvc_k2",  "2 2 1 0 0 0",(2,2,1),150,50,400),
 ("conv_yvc_k4",  "4 4 2 0 0 0",(4,4,2),150,50,400),   # baseline
 ("conv_yvc_nqx1","4 4 2 0 0 0",(1,1,1),150,50,400),
 ("conv_yvc_nqx2","4 4 2 0 0 0",(2,2,1),150,50,400),
 ("conv_yvc_ef100","4 4 2 0 0 0",(4,4,2),100,50,400),
 ("conv_yvc_ef200","4 4 2 0 0 0",(4,4,2),200,50,400),
 ("conv_yvc_ec40","4 4 2 0 0 0",(4,4,2),150,40,320),
 ("conv_yvc_ec60","4 4 2 0 0 0",(4,4,2),150,60,480),
]
man=[]
for tag,k,nqx,ef,ew,er in SWEEP:
    man.append(emit(tag,cell,atoms,k,nqx,ef,ew,er,0.02))
    print(f"  {tag:14s} k={k.split()[0:3]} nqx={nqx} ecutfock={ef} ew/er={ew}/{er}")
open(f"{OUT}/conv_manifest.txt","w").write("\n".join(man)+"\n")
print(f"\nwrote {len(man)} convergence-sweep inputs (Y-vc loaded)")
