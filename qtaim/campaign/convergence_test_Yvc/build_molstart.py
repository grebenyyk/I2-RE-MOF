#!/usr/bin/env python3
"""Construct molecular-START geometries (4 I2, inter-pair >3.5A) for vc/st by
re-partnering the polyiodide I into 4 pairs and spacing the pair-centers ~6.3A apart
along the polyiodide principal axis. Forces the molecular basin on relax.
Framework frozen, I free. Emits PBE/4x4x2 relax inputs mmol_Y{v}.in."""
import numpy as np
OUT="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
PSE={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF","Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
MASS={"C":12.011,"H":1.008,"O":15.999,"Y":88.906,"I":126.904}
def load(tag):
    L=open(f"{OUT}/{tag}.in").read().splitlines()
    ci=[i for i,l in enumerate(L) if l.startswith("CELL_PARAMETERS")][0]
    ai=[i for i,l in enumerate(L) if l.startswith("ATOMIC_POSITIONS")][0]
    cell=np.array([[float(x) for x in L[ci+1+k].split()[:3]] for k in range(3)])
    fw=[];If=[];k=ai+1
    while k<len(L):
        s=L[k].split()
        if not s or s[0] in("K_POINTS","ATOMIC_SPECIES"): break
        if len(s)>=4 and s[0] in PSE:
            f=[float(s[1]),float(s[2]),float(s[3])]
            (If if s[0]=="I" else fw).append((s[0],f))
        k+=1
    return cell,fw,np.array([x for _,x in If])
def greedy_pair(P):
    used=set();pairs=[];n=len(P)
    ds=sorted([(np.linalg.norm(P[i]-P[j]),i,j) for i in range(n) for j in range(i+1,n)])
    for d,i,j in ds:
        if i not in used and j not in used: pairs.append((i,j));used|={i,j}
        if len(pairs)==4: break
    return pairs
def build(v,spacing=6.3):
    cell,fw,Ifrac=load(f"rd_poly_Y{v}")
    I=Ifrac@cell
    pairs=greedy_pair(I)
    # set each pair to 2.70 A around its midpoint
    for (i,j) in pairs:
        c=(I[i]+I[j])/2; u=I[j]-I[i]; u=u/np.linalg.norm(u)
        I[i]=c-1.35*u; I[j]=c+1.35*u
    # principal axis (channel direction)
    cen=I.mean(0);_,_,Vt=np.linalg.svd(I-cen);ax=Vt[0]
    pc=np.array([(I[i]+I[j])/2 for (i,j) in pairs])
    order=np.argsort(np.dot(pc-cen,ax))
    new=np.array([-1.5,-0.5,0.5,1.5])*spacing  # 4 centers spaced 'spacing' apart, centered
    for k,idx in enumerate(order):
        (i,j)=pairs[idx];c=(I[i]+I[j])/2
        shift=(new[k]-np.dot(c-cen,ax))*ax
        I[i]+=shift;I[j]+=shift
    Ifrac2=(I@np.linalg.inv(cell))%1.0
    # checks
    ii=[]
    for i in range(8):
        for j in range(i+1,8):
            df=I[i]-I[j];df-=np.round(df);d=np.linalg.norm(df@cell)
            if d<3.6: ii.append(round(d,2))
    fwc=np.array([f@cell for _,f in fw])
    minIF=min(np.linalg.norm((fwc-p)-np.round(fwc-p)@cell,axis=1).min() for p in I)
    nshort=len([d for d in ii if d<2.9])
    print(f"{v}: pairs(sites)={[(i+1,j+1) for i,j in pairs]}")
    print(f"  I-I<3.6: {sorted(ii)}  -> {nshort} pair(s) <2.9A ({'OK 4 I2' if nshort==4 else 'CHECK'})")
    print(f"  min I-framework = {minIF:.2f} A  ({'no clash' if minIF>2.4 else 'CLASH!'})")
    # emit PBE relax input (framework frozen if_pos 000, I free 111)
    atoms=[(e,f) for e,f in fw]+[("I",list(Ifrac2[k])) for k in range(8)]
    uniq=[];[uniq.append(a[0]) for a in atoms if a[0] not in uniq]
    nfw=len(fw)
    L=[f"&CONTROL\n  calculation='relax'\n  prefix='mmol_Y{v}'\n  outdir='./tmp_mmol_Y{v}'",
       "  pseudo_dir='/workspace/qe-7.5/pseudo'\n  tprnfor=.true.\n  etot_conv_thr=1.0d-5\n  forc_conv_thr=5.0d-3\n/",
       f"&SYSTEM\n  ibrav=0\n  nat={len(atoms)}\n  ntyp={len(uniq)}\n  ecutwfc=50\n  ecutrho=400",
       "  occupations='smearing'\n  smearing='mv'\n  degauss=0.02\n  nspin=1\n  vdw_corr='grimme-d3'\n  dftd3_version=4\n/",
       "&ELECTRONS\n  conv_thr=1.0d-7\n  diagonalization='david'\n  diago_david_ndim=2",
       "  mixing_mode='local-TF'\n  mixing_beta=0.2\n  mixing_ndim=24\n  electron_maxstep=300\n/",
       "&IONS\n  ion_dynamics='bfgs'\n  ion_maxstep=80\n/","ATOMIC_SPECIES"]
    for e in uniq: L.append(f"  {e:<2s}{MASS[e]:10.3f}  {PSE[e]}")
    L.append("CELL_PARAMETERS angstrom")
    for r in cell: L.append(f"  {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
    L.append("ATOMIC_POSITIONS crystal")
    for k,(e,f) in enumerate(atoms):
        fl="0 0 0" if k<nfw else "1 1 1"
        L.append(f"  {e:<2s} {f[0]:.10f} {f[1]:.10f} {f[2]:.10f}  {fl}")
    L.append("K_POINTS automatic\n  4 4 2 0 0 0\n")
    open(f"{OUT}/mmol_Y{v}.in","w").write("\n".join(L))
    print(f"  wrote mmol_Y{v}.in (forc_conv_thr=5e-3, ion_maxstep=80)")
for v in ["vc","st"]: build(v)
