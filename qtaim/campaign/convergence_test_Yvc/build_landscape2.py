#!/usr/bin/env python3
"""Landscape enrichment (batch 2): 3 more void-placed molecular (S4,S5,S6) +
3 more perturbed polyiodide (S4,S5,S6) per variant. S6 molecular uses a-axis I2
orientation (vs channel/c in S1-5) for orientational diversity. PBE/4x4x2 relax.
st molecular -> mixing_beta 0.05 (known SCF-stall fix); else 0.2."""
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
def emit_relax(tag,cell,fw,Ifrac,beta):
    aa=list(fw)+[("I",*p) for p in Ifrac]; uniq=[];[uniq.append(a[0]) for a in aa if a[0] not in uniq]; nfw=len(fw)
    LL=[f"&CONTROL\n  calculation='relax'\n  prefix='{tag}'\n  outdir='./tmp_{tag}'","  pseudo_dir='/workspace/qe-7.5/pseudo'\n  tprnfor=.true.\n  etot_conv_thr=1.0d-5\n  forc_conv_thr=5.0d-3\n/",
       f"&SYSTEM\n  ibrav=0\n  nat={len(aa)}\n  ntyp={len(uniq)}\n  ecutwfc=50\n  ecutrho=400","  occupations='smearing'\n  smearing='mv'\n  degauss=0.02\n  nspin=1\n  vdw_corr='grimme-d3'\n  dftd3_version=4\n/",
       f"&ELECTRONS\n  conv_thr=1.0d-7\n  diagonalization='david'\n  diago_david_ndim=2\n  mixing_mode='local-TF'\n  mixing_beta={beta}\n  mixing_ndim=24\n  electron_maxstep=300\n/","&IONS\n  ion_dynamics='bfgs'\n  ion_maxstep=80\n/","ATOMIC_SPECIES"]
    for e in uniq: LL.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    LL+=["CELL_PARAMETERS angstrom"]+[f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}" for r in cell]+["ATOMIC_POSITIONS crystal"]
    for k,a in enumerate(aa): LL.append(f"  {a[0]:<2s} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}  {'0 0 0' if k<nfw else '1 1 1'}")
    LL.append("K_POINTS automatic\n  4 4 2 0 0 0\n")
    open(f"{OUT}/{tag}.in","w").write("\n".join(LL))
# ---- molecular: diversify via perturbation of TWO existing minima (S4,S5) +
#      one void-placement with a-axis orientation (S6). Robust (valid starts) + diverse. ----
MOL_SRC_BEST={"act":"rd_mol_Yact","vc":"hmol_Yvc","st":"hse_lmol_Yst_S3"}   # lowest-E molecular per variant
MOL_SRC_ALT ={"act":"hse_lmol_Yact_S2","vc":"hse_lmol_Yvc_S2","st":"hse_lmol_Yst_S2"}  # higher-E molecular
def load_I_from(srctag):
    L=open(f"{OUT}/{srctag}.in").read().splitlines()
    ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][0]; I=[];k=ai+1
    while k<len(L):
        s=L[k].split()
        if not s or s[0] in("K_POINTS","ATOMIC_SPECIES"): break
        if len(s)>=4 and s[0]=="I": I.append([float(s[1]),float(s[2]),float(s[3])])
        k+=1
    return np.array(I)
def void_place(cell,inv,fwc,axis):
    aL,bL,cL=[np.linalg.norm(cell[i]) for i in range(3)]
    # rank c-layers by post-placement clearance
    rank=[]
    for cf in np.arange(0.10,0.91,0.05):
        best=-1
        for af in np.arange(0.05,0.96,0.05):
            for bf in np.arange(0.05,0.96,0.05):
                if axis=='c': pts=[(af,bf,(cf*cL-1.35)/cL),(af,bf,(cf*cL+1.35)/cL)]
                else: pts=[((af*aL-1.35)/aL,bf,cf),((af*aL+1.35)/aL,bf,cf)]
                cc=np.array([p@cell for p in pts]); dm=min(np.linalg.norm(((fwc-p)@inv-np.round((fwc-p)@inv))@cell,axis=1).min() for p in cc)
                if dm>best: best=dm
        if best>2.4: rank.append((best,cf))
    rank.sort(reverse=True)
    # two best-separated layers
    chosen=[]
    for d,cf in rank:
        if all(abs(cf-q)>0.3 for q in chosen): chosen.append(cf)
        if len(chosen)>=2: break
    if len(chosen)<2: return None
    I=[]; minIF=9
    for cf in chosen:
        cand=[]
        for af in np.arange(0.05,0.96,0.05):
            for bf in np.arange(0.05,0.96,0.05):
                if axis=='c': pts=[(af,bf,(cf*cL-1.35)/cL),(af,bf,(cf*cL+1.35)/cL)]
                else: pts=[((af*aL-1.35)/aL,bf,cf),((af*aL+1.35)/aL,bf,cf)]
                cc=np.array([p@cell for p in pts]); dm=min(np.linalg.norm(((fwc-p)@inv-np.round((fwc-p)@inv))@cell,axis=1).min() for p in cc)
                if dm>2.4: cand.append((dm,af,bf))
        cand.sort(reverse=True)
        for d,af,bf in cand[:2]:
            if axis=='c': I+=([[af,bf,(cf*cL-1.35)/cL],[af,bf,(cf*cL+1.35)/cL]])
            else: I+=([[(af*aL-1.35)/aL,bf,cf],[(af*aL+1.35)/aL,bf,cf]])
            minIF=min(minIF,d)
    return (np.array(I)[:8],chosen,minIF) if len(I)>=8 else None
for v in ["act","vc","st"]:
    cell,fw=load_fw(f"rd_poly_Y{v}"); inv=np.linalg.inv(cell)
    fwc=np.array([[x,y,z]@cell for _,x,y,z in fw]); beta=0.05 if v=="st" else 0.2
    # S4: perturb BEST molecular minimum
    I=load_I_from(MOL_SRC_BEST[v]); rng=np.random.default_rng(401)
    Ip=(I+rng.uniform(-0.03,0.03,I.shape))%1.0; emit_relax(f"lmol_Y{v}_S4",cell,fw,Ip,beta)
    print(f"[mol] lmol_Y{v}_S4: perturb {MOL_SRC_BEST[v]} (±0.03), beta={beta}")
    # S5: perturb ALT molecular minimum (different basin point)
    I=load_I_from(MOL_SRC_ALT[v]); rng=np.random.default_rng(402)
    Ip=(I+rng.uniform(-0.03,0.03,I.shape))%1.0; emit_relax(f"lmol_Y{v}_S5",cell,fw,Ip,beta)
    print(f"[mol] lmol_Y{v}_S5: perturb {MOL_SRC_ALT[v]} (±0.03), beta={beta}")
    # S6: void-placement, a-axis orientation (fallback c if clashes)
    res=void_place(cell,inv,fwc,'a')
    if res is None or res[2]<2.4:
        res=void_place(cell,inv,fwc,'c'); ax='c(fallback)'
    else: ax='a'
    Ip,clays,minIF=res; emit_relax(f"lmol_Y{v}_S6",cell,fw,Ip,beta)
    print(f"[mol] lmol_Y{v}_S6: void-place axis={ax}, c-layers={[round(c,2) for c in clays]}, min I-FW={minIF:.2f}A, beta={beta}")
# ---- polyiodide: (variant, seed, rng_seed, perturbation) ----
POLY=[("act",4,301,0.03),("act",5,302,0.04),("act",6,303,0.02),
      ("vc",4,311,0.03),("vc",5,312,0.04),("vc",6,313,0.02),
      ("st",4,321,0.03),("st",5,322,0.04),("st",6,323,0.02)]
for v,seed,rs,kick in POLY:
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
        a2[j][1]=(a2[j][1]+rng.uniform(-kick,kick))%1.0
        a2[j][2]=(a2[j][2]+rng.uniform(-kick,kick))%1.0
        a2[j][3]=(a2[j][3]+rng.uniform(-kick,kick))%1.0
    tag=f"lpoly_Y{v}_S{seed}"; emit_relax(tag,cell,fw,np.array([a[1:] for a in a2[nfw:]]),0.2)
    print(f"[poly] {tag}: kick=±{kick}, rng_seed={rs}")
print("\nbatch-2 manifest:")
for v in ["act","vc","st"]:
    for k in ["lmol","lpoly"]:
        tags=sorted(f for f in os.listdir(OUT) if f.startswith(f"{k}_Y{v}_S") and f.endswith(".in") and int(f.split("_S")[1].split('.')[0])>=4)
        print(f"  {k}_Y{v} (new): {[t.replace('.in','') for t in tags]}")
