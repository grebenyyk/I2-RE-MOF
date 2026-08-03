#!/usr/bin/env python3
"""Fast vectorized Hirshfeld population analysis on a QE valence cube.
Density-based partition (no projection basis -> independent of Mulliken/projwfc).
Uses free-atom valence densities from each PP's PP_RHOATOM (P(r)=4pir^2 rho).
Reports per-atom Hirshfeld population and charge q = Z_val - Pop."""
import numpy as np, re, sys, os
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
PSE="/workspace/qe-7.5/pseudo"
PPMAP={"C":"C.pbe-n-kjpaw_psl.1.0.0.UPF","H":"H.pbe-kjpaw.UPF","O":"O.pbe-kjpaw.UPF",
       "Y":"Y.pbe-spn-kjpaw_psl.1.0.0.UPF","Tb":"Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF",
       "Gd":"Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF","I":"I.pbe-n-kjpaw_psl.1.0.0.UPF"}
ZVAL={"C":4,"H":1,"O":6,"Y":11,"Tb":11,"Gd":11,"I":7}
PI4=4*np.pi
_ppcache={}
def free_atom(species):
    if species in _ppcache: return _ppcache[species]
    t=open(f"{PSE}/{PPMAP[species]}").read()
    # PP_R mesh
    rm=re.search(r"<PP_R[^>]*>(.*?)</PP_R>",t,re.S)
    r=np.array([float(x) for x in rm.group(1).split()])
    rh=re.search(r"<PP_RHOATOM[^>]*>(.*?)</PP_RHOATOM>",t,re.S)
    P=np.array([float(x) for x in rh.group(1).split()])  # P(r)=4pi r^2 rho(r), integral=Zval
    rho=P/(PI4*r*r+1e-12); rho[0]=P[1]/(PI4*r[1]*r[1]) if r[0]==0 else rho[0]
    _trap=getattr(np,'trapezoid',None) or getattr(np,'trapz')
    chk=_trap(P,r)
    print(f"  [pp] {species}: {len(r)} mesh pts, integral(P)={chk:.4f} (Zval={ZVAL[species]})",flush=True)
    _ppcache[species]=(r,rho)
    return r,rho
def read_cube(fn):
    L=open(fn).readlines(); nat=int(L[2].split()[0])
    nx,ny,nz=[int(L[i].split()[0]) for i in (3,4,5)]
    ga=np.array([float(x) for x in L[3].split()[1:4]]); gb=np.array([float(x) for x in L[4].split()[1:4]]); gc=np.array([float(x) for x in L[5].split()[1:4]])
    d=np.array([float(x) for l in L[6+nat:] for x in l.split()]).reshape((nx,ny,nz))
    return d,nx,ny,nz,np.array([ga*nx,gb*ny,gc*nz]).T  # cellB(bohr) cols
def read_atoms_scf(scf):
    s=open(scf).read()
    cellA=np.array([list(map(float,re.findall(r"-?\d+\.\d+",l))) for l in s.split("CELL_PARAMETERS")[1].splitlines()[1:4]])
    apos=[]
    for l in s.split("ATOMIC_POSITIONS")[1].split("K_POINTS")[0].splitlines():
        m=re.match(r"\s*([A-Z][a-z]?)\s+(-?\d.\d+)\s+(-?\d.\d+)\s+(-?\d.\d+)",l)
        if m: apos.append((m.group(1),np.array([float(m.group(2)),float(m.group(3)),float(m.group(4))])))
    return cellA,apos
def hirshfeld_cube(cube,scf):
    d,nx,ny,nz,cellB=read_cube(cube)
    cellA,apos=read_atoms_scf(scf)
    Vcell=np.linalg.det(cellB); Vvox=Vcell/(nx*ny*nz)
    ii,jj,kk=np.meshgrid(np.arange(nx),np.arange(ny),np.arange(nz),indexing='ij')
    fx=ii/nx; fy=jj/ny; fz=kk/nz
    rho=d
    promol=np.zeros_like(rho)
    for sp,fc in apos:
        r,rhoA=free_atom(sp)
        dfx=fx-fc[0]; dfx=dfx-np.round(dfx); dfy=fy-fc[1]; dfy=dfy-np.round(dfy); dfz=fz-fc[2]; dfz=dfz-np.round(dfz)
        cx=dfx*cellB[0,0]+dfy*cellB[1,0]+dfz*cellB[2,0]; cy=dfx*cellB[0,1]+dfy*cellB[1,1]+dfz*cellB[2,1]; cz=dfx*cellB[0,2]+dfy*cellB[1,2]+dfz*cellB[2,2]
        dist=np.clip(np.sqrt(cx*cx+cy*cy+cz*cz),r[0],r[-1])
        promol+=np.interp(dist,r,rhoA)
    promol=np.maximum(promol,1e-12)
    out=[]
    for sp,fc in apos:
        r,rhoA=free_atom(sp)
        dfx=fx-fc[0]; dfx=dfx-np.round(dfx); dfy=fy-fc[1]; dfy=dfy-np.round(dfy); dfz=fz-fc[2]; dfz=dfz-np.round(dfz)
        cx=dfx*cellB[0,0]+dfy*cellB[1,0]+dfz*cellB[2,0]; cy=dfx*cellB[0,1]+dfy*cellB[1,1]+dfz*cellB[2,1]; cz=dfx*cellB[0,2]+dfy*cellB[1,2]+dfz*cellB[2,2]
        dist=np.clip(np.sqrt(cx*cx+cy*cy+cz*cz),r[0],r[-1])
        w=np.interp(dist,r,rhoA)/promol
        out.append((sp,fc,float((rho*w).sum()*Vvox)))
    return out
def diff_hirshfeld(sys):
    """delta = Pop(loaded)-Pop(empty)-Pop(i2only), atoms matched by INDEX.
    loaded = Cx32,Hx12,Ox24,Yx4,Ix8 (80). empty = framework only (first 72, same order).
    i2only = Ix8 (last 8 of loaded, same order)."""
    L=hirshfeld_cube(f"{BASE}/{sys}/val.cube",f"{BASE}/{sys}/scf_loaded.in")
    E=hirshfeld_cube(f"{BASE}/{sys}/empty_val.cube",f"{BASE}/{sys}/scf_empty.in")
    I=hirshfeld_cube(f"{BASE}/{sys}/i2only_val.cube",f"{BASE}/{sys}/scf_i2only.in")
    nE=len(E); nI=len(I)
    rows=[]
    for i,(sp,fc,pL) in enumerate(L):
        pE = E[i][2] if i<nE else 0.0          # framework atom matches empty by index
        pI = I[i-nE][2] if i>=nE and (i-nE)<nI else 0.0  # I atom matches i2only by index
        rows.append((sp,pL,pE,pI,pL-pE-pI))
    return rows
if __name__=="__main__":
    from collections import defaultdict
    for sys in (sys.argv[1:] or ["Y_str","Gd_str"]):
        print(f"=== differential Hirshfeld {sys}: deltaPop = Pop(loaded)-Pop(empty)-Pop(i2only) ===",flush=True)
        rows=diff_hirshfeld(sys)
        agg=defaultdict(lambda:[0.0,0])
        print(f"  {'sp':<3}{'Pop_L':>8}{'Pop_E':>8}{'Pop_I2':>8}{'delta':>9}")
        for sp,pL,pE,pI,dp in rows:
            agg[sp][0]+=dp; agg[sp][1]+=1
            if sp=="I": print(f"  {sp:<3}{pL:>8.4f}{pE:>8.4f}{pI:>8.4f}{dp:>+9.4f}")
        print("  --- sum deltaPop by species (e/cell; + = atom GAINED = reduced) ---")
        tot=0
        for sp in ["C","H","O","Y","Tb","Gd","I"]:
            if agg[sp][1]: print(f"  {sp:<3} N={agg[sp][1]:>3}  sum deltaPop = {agg[sp][0]:>+.4f} e  (per-atom {agg[sp][0]/agg[sp][1]:>+.4f})"); tot+=agg[sp][0]
        print(f"  TOTAL deltaPop (sum check, must ~0): {tot:+.4f} e")

