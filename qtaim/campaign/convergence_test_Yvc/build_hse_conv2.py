#!/usr/bin/env python3
"""Round 2: ecutwfc binding check + 6x6x3 k-check on Y-vc.
ecutwfc check: Y-vc empty + gas at 40/320, 60/480 (50/400 already done:
  loaded=conv_yvc_k4, empty=hse_empty_yvc, gas=hse_i2_vac) -> binding trend.
6x6x3: Y-vc loaded at k=6x6x3, nqx=6x6x3, ef100, 50/400 vs the 4x4x2 baseline."""
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
def emit(tag,cell,atoms,kmesh,nqx,ef,ew,er,dg,smode="mv"):
    uniq=[]; [uniq.append(e) for e in [a[0] for a in atoms] if e not in uniq]
    L=[f"&CONTROL\n  calculation = 'scf'\n  prefix = '{tag}'\n  outdir = './tmp_{tag}'",
       "  pseudo_dir = '/workspace/qe-7.5/pseudo'\n  tprnfor = .true.\n/",
       f"&SYSTEM\n  ibrav = 0\n  nat = {len(atoms)}\n  ntyp = {len(uniq)}\n  ecutwfc = {ew}\n  ecutrho = {er}",
       f"  occupations = 'smearing'\n  smearing = '{smode}'\n  degauss = {dg}\n  nspin = 1",
       "  input_dft = 'hse'\n  exxdiv_treatment = 'gygi-baldereschi'\n  x_gamma_extrapolation = .true.",
       f"  exx_fraction = 0.25\n  screening_parameter = 0.106\n  ecutfock = {ef}",
       f"  nqx1 = {nqx[0]}\n  nqx2 = {nqx[1]}\n  nqx3 = {nqx[2]}",
       "  vdw_corr = 'grimme-d3'\n  dftd3_version = 4\n/",
       "&ELECTRONS\n  conv_thr = 1.0d-6\n  diagonalization = 'david'\n  diago_david_ndim = 2",
       f"  mixing_mode = 'local-TF'\n  mixing_beta = 0.05\n  mixing_ndim = 24\n  electron_maxstep = 300\n/",
       "ATOMIC_SPECIES"]
    for e in uniq: L.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell: L.append(f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for e,x,y,z in atoms: L.append(f"  {e:<2s} {x:.10f} {y:.10f} {z:.10f}")
    L.append(f"K_POINTS automatic\n  {kmesh}\n")
    open(f"{OUT}/{tag}.in","w").write("\n".join(L)); return tag
ce,ae=cif_geom(f"{CIF}/vc-relaxed/Y_vc_empty.cif")          # empty (insulator -> Gamma)
cg,ag=inp_geom(f"{OUT}/hse_i2_vac.in")                       # gas (molecule -> Gamma)
cl,al=cif_geom(f"{CIF}/vc-relaxed/Y_vc_loaded.cif")          # loaded for 6x6x3
man=[]
# ecutwfc binding check: empty + gas at 40/60 (50 already have), Gamma, nqx111
for ew,er,tt in [(40,320,"ec40"),(60,480,"ec60")]:
    man.append(emit(f"chk_empty_yvc_{tt}",ce,ae,"1 1 1 0 0 0",(1,1,1),100,ew,er,0.01))
    man.append(emit(f"chk_gas_{tt}",cg,ag,"1 1 1 0 0 0",(1,1,1),100,ew,er,0.01))
# 6x6x3 k-check: loaded, nqx=6x6x3, ef100, 50/400
man.append(emit("chk_yvc_k663",cl,al,"6 6 3 0 0 0",(6,6,3),100,50,400,0.02))
for t in man: print("wrote",t)
open(f"{OUT}/conv2_manifest.txt","w").write("\n".join(man)+"\n")
