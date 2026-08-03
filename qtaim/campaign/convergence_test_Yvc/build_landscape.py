#!/usr/bin/env python3
"""Expanded landscape: 2 new void-placed molecular (S2,S3, different c-layers)
+ 2-3 new perturbed polyiodide per variant. Existing: mol S1 (act=rd_mol,vc=hmol,st=mmol),
poly S1 (vc,st=hmol_S1). All PBE/4x4x2 relax -> HSE/4x4x2 SP."""
import numpy as np, re, os
OUT="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF","Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"I":126.904}
def load_fw(tag):
    L=open(f"{OUT}/{tag}.in").read().splitlines()
    ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
    ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][0]
    cell=np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)]); fw=[];k=ai+1
    while k<len(L):
        s=L[k].split()
        if not s or s[0] in("K_POINTS","ATOMIC_SPECIES"): break
        if len(s)>=4 and s[0] in("C","H","O","Y"): fw.append((s[0],float(s[1]),float(s[2]),float(s[3])))
        k+=1
    return cell,fw
def mindist(inv,cell,fwc,p):
    fr=(fwc-p)@inv; fr-=np.round(fr); return np.linalg.norm(fr@cell,axis=1).min()
def emit_relax(tag,cell,fw,Ifrac):
    aa=list(fw)+[("I",*p) for p in Ifrac]; uniq=[];[uniq.append(a[0]) for a in aa if a[0] not in uniq]; nfw=len(fw)
    LL=[f"&CONTROL\n  calculation='relax'\n  prefix='{tag}'\n  outdir='./tmp_{tag}'","  pseudo_dir='/workspace/qe-7.5/pseudo'\n  tprnfor=.true.\n  etot_conv_thr=1.0d-5\n  forc_conv_thr=5.0d-3\n/",
       f"&SYSTEM\n  ibrav=0\n  nat={len(aa)}\n  ntyp={len(uniq)}\n  ecutwfc=50\n  ecutrho=400","  occupations='smearing'\n  smearing='mv'\n  degauss=0.02\n  nspin=1\n  vdw_corr='grimme-d3'\n  dftd3_version=4\n/",
       "&ELECTRONS\n  conv_thr=1.0d-7\n  diagonalization='david'\n  diago_david_ndim=2\n  mixing_mode='local-TF'\n  mixing_beta=0.2\n  mixing_ndim=24\n  electron_maxstep=300\n/","&IONS\n  ion_dynamics='bfgs'\n  ion_maxstep=80\n/","ATOMIC_SPECIES"]
    for e in uniq: LL.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    LL+=["CELL_PARAMETERS angstrom"]+[f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}" for r in cell]+["ATOMIC_POSITIONS crystal"]
    for k,a in enumerate(aa): LL.append(f"  {a[0]:<2s} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}  {'0 0 0' if k<nfw else '1 1 1'}")
    LL.append("K_POINTS automatic\n  4 4 2 0 0 0\n")
    open(f"{OUT}/{tag}.in","w").write("\n".join(LL))
# --- Molecular (void-placed, varying c-layers) ---
MOL_CLAYS=[("act",2,[0.15,0.65]),("act",3,[0.35,0.85]),("vc",2,[0.15,0.65]),("vc",3,[0.35,0.85]),("st",2,[0.15,0.65]),("st",3,[0.35,0.85])]
for v,seed,clays in MOL_CLAYS:
    cell,fw=load_fw(f"rd_poly_Y{v}"); inv=np.linalg.inv(cell); fwc=np.array([f@cell for _,*f in [(e,x,y,z) for e,x,y,z in fw]]); fwc=np.array([[x,y,z]@cell for _,x,y,z in fw])
    aL,bL,cL=[np.linalg.norm(cell[i]) for i in range(3)]; Iatoms=[]
    for cf in clays:
        cand=[]
        for af in np.arange(0.05,0.96,0.05):
            for bf in np.arange(0.05,0.96,0.05):
                d=mindist(inv,cell,fwc,np.array([af,bf,cf])@cell)
                if d>2.5: cand.append((d,af,bf))
        cand.sort(reverse=True); s0=cand[0] if cand else None; s1=None
        if s0:
            p0=np.array([s0[1],s0[2]])@cell[:2,:2]
            for d,a,b in cand[1:]:
                if np.linalg.norm(np.array([a,b])@cell[:2,:2]-p0)>4.0: s1=(d,a,b); break
        for d,af,bf in [s0,s1] if s0 and s1 else []:
            Iatoms.append([af,bf,(cf*cL-1.35)/cL]); Iatoms.append([af,bf,(cf*cL+1.35)/cL])
    Ifrac=np.array(Iatoms)%1.0; tag=f"lmol_Y{v}_S{seed}"; emit_relax(tag,cell,fw,Ifrac)
    ni=sum(1 for p in Ifrac); print(f"[mol] {tag}: {ni} I, c-layers={clays}")
# --- Polyiodide (perturbed, varying seeds) ---
POLY=[("act",1,201),("act",2,202),("act",3,203),("vc",2,212),("vc",3,213),("st",2,222),("st",3,223)]
for v,seed,rs in POLY:
    cell,fw=load_fw(f"rd_poly_Y{v}")
    L=open(f"{OUT}/rd_poly_Y{v}.in").read().splitlines()
    ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][0]; atoms=[];k=ai+1
    while k<len(L):
        s=L[k].split()
        if not s or s[0] in("K_POINTS","ATOMIC_SPECIES"): break
        if len(s)>=4 and s[0] in PSE:
            try: atoms.append((s[0],float(s[1]),float(s[2]),float(s[3])))
            except: pass
        k+=1
    nfw=len(atoms)-8; rng=np.random.default_rng(rs); a2=[list(a) for a in atoms]
    for j in range(nfw,len(atoms)):
        a2[j][1]=(a2[j][1]+rng.uniform(-0.02,0.02))%1.0; a2[j][2]=(a2[j][2]+rng.uniform(-0.02,0.02))%1.0; a2[j][3]=(a2[j][3]+rng.uniform(-0.02,0.02))%1.0
    tag=f"lpoly_Y{v}_S{seed}"; emit_relax(tag,cell,fw,np.array([a[1:] for a in a2[nfw:]]))
    print(f"[poly] {tag}: perturbed polyiodide, kick_seed={rs}")
print(f"\nmanifest:")
for v in ["act","vc","st"]:
    for k in ["lmol","lpoly"]:
        tags=[f for f in os.listdir(OUT) if f.startswith(f"{k}_Y{v}_") and f.endswith(".in")]
        print(f"  {k}_Y{v}: {sorted(tags)}")
