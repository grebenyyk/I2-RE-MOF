#!/usr/bin/env python3
"""Uniform 9-system HSE06 redo. ALL geometries from authoritative CIFs.
Locked settings (from Y-vc sweep): loaded 4x4x2 (metallic, converged); empty+gas Gamma
(insulator/molecule, k-exact); degauss 0.02 mv (UNIFORM); ecutfock 100; ecutwfc/rho 50/400;
nqx default(=k); HSE06+D3(BJ)v4; f-in-core PPs, nspin=1.
conv_thr 1e-6 loaded/gas, 1e-3 empty (vc-empty limit-cycle guard, energy converged <1kJ).
Replaces the bug-afflicted table (#42): CIF geometries fix Gd-act + vc-empty; uniform mesh
fixes mixed Gamma/4x4x2 + Gamma-only-Fock on act; fresh consistent gas fixes the -13kJ ref."""
import re,math,numpy as np,os
BASE="/workspace/dimitry/calc/qe/i2_mof"; CIF=f"{BASE}/final_activated_results"; OUT=f"{BASE}/qtaim/campaign/convergence_test_Yvc"
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
     "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF","Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"Gd":157.25,"Tb":158.925,"I":126.904}
def cif_geom(f):
    t=open(f).read()
    def g(p): return float(re.search(p,t).group(1))
    a,b,c=[g(rf"_cell_length_{n}\s+([\d.\-eE]+)") for n in("a","b","c")]
    al,be,ga=[math.radians(g(rf"_cell_angle_{n}\s+([\d.\-eE]+)")) for n in("alpha","beta","gamma")]
    cell=np.zeros((3,3)); cell[0]=[a,0,0]; cell[1]=[b*math.cos(ga),b*math.sin(ga),0]
    cg=(math.cos(al)-math.cos(be)*math.cos(ga))/math.sin(ga)
    cell[2]=[c*math.cos(be),c*cg,c*math.sqrt(max(0,1-math.cos(be)**2-cg**2))]
    atoms=[]
    for l in t.splitlines():
        s=l.split()
        if len(s)>=5 and s[1] in PSE:
            try: atoms.append((s[1],float(s[2]),float(s[3]),float(s[4])))
            except: pass
    return cell,atoms
def inp_geom(f):
    L=open(f).read().splitlines()
    ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
    ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][0]
    cell=np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)]); atoms=[];k=ai+1
    while k<len(L):
        s=L[k].split()
        if not s or s[0] in("K_POINTS","ATOMIC_SPECIES"): break
        if len(s)>=4 and s[0] in PSE:
            try: atoms.append((s[0],float(s[1]),float(s[2]),float(s[3])))
            except: pass
        k+=1
    return cell,atoms
def emit(tag,cell,atoms,kmesh,nqx,ef,ew,er,dg,convthr):
    uniq=[]; [uniq.append(e) for e in [a[0] for a in atoms] if e not in uniq]
    L=[f"&CONTROL\n  calculation = 'scf'\n  prefix = '{tag}'\n  outdir = './tmp_{tag}'",
       "  pseudo_dir = '/workspace/qe-7.5/pseudo'\n  tprnfor = .true.\n/",
       f"&SYSTEM\n  ibrav = 0\n  nat = {len(atoms)}\n  ntyp = {len(uniq)}\n  ecutwfc = {ew}\n  ecutrho = {er}",
       f"  occupations = 'smearing'\n  smearing = 'mv'\n  degauss = {dg}\n  nspin = 1",
       "  input_dft = 'hse'\n  exxdiv_treatment = 'gygi-baldereschi'\n  x_gamma_extrapolation = .true.",
       f"  exx_fraction = 0.25\n  screening_parameter = 0.106\n  ecutfock = {ef}",
       f"  nqx1 = {nqx[0]}\n  nqx2 = {nqx[1]}\n  nqx3 = {nqx[2]}",
       "  vdw_corr = 'grimme-d3'\n  dftd3_version = 4\n/",
       f"&ELECTRONS\n  conv_thr = {convthr}\n  diagonalization = 'david'\n  diago_david_ndim = 2",
       "  mixing_mode = 'local-TF'\n  mixing_beta = 0.05\n  mixing_ndim = 24\n  electron_maxstep = 300\n/",
       "ATOMIC_SPECIES"]
    for e in uniq: L.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell: L.append(f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for e,x,y,z in atoms: L.append(f"  {e:<2s} {x:.10f} {y:.10f} {z:.10f}")
    L.append(f"K_POINTS automatic\n  {kmesh}\n")
    open(f"{OUT}/{tag}.in","w").write("\n".join(L)); return tag
SUB={"st":"solvent-induced","act":"activated","vc":"vc-relaxed"}
METALS=[("Y","Y"),("Gd","Gd"),("Tb","Tb")]
man=[]
for ml,uc in METALS:
    for v in ["st","act","vc"]:
        cl,al=cif_geom(f"{CIF}/{SUB[v]}/{uc}_{v}_loaded.cif")
        man.append(emit(f"rd_poly_{ml}{v}",cl,al,"4 4 2 0 0 0",(4,4,2),100,50,400,0.02,"1.0d-6"))
        ce,ae=cif_geom(f"{CIF}/{SUB[v]}/{uc}_{v}_empty.cif")
        man.append(emit(f"rd_empty_{ml}{v}",ce,ae,"1 1 1 0 0 0",(1,1,1),100,50,400,0.02,"1.0d-3"))
cg,ag=inp_geom(f"{OUT}/chk_gas_ec50.in")
man.append(emit("rd_gas",cg,ag,"1 1 1 0 0 0",(1,1,1),100,50,400,0.02,"1.0d-6"))
for t in man:
    nat=len([l for l in open(f"{OUT}/{t}.in") if l.startswith("  I ") or l.startswith("  C ") or l.startswith("  H ") or l.startswith("  O ") or l.startswith("  Y ") or l.startswith("  Gd ") or l.startswith("  Tb ")])
    nI=len([l for l in open(f"{OUT}/{t}.in") if l.startswith("  I ")])
    print(f"  {t:16s} atoms~{nat:3d} I={nI}")
open(f"{OUT}/redo_manifest.txt","w").write("\n".join(man)+"\n")
print(f"\nwrote {len(man)} redo inputs (9 loaded 4x4x2 + 9 empty Gamma + 1 gas)")
