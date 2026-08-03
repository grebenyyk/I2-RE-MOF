#!/usr/bin/env python3
"""Per-system real-space analysis: BCP+EML (Abramov), CT table (Bader basins), ELF slice.
Run after SCFs. Usage: python3 analyze_campaign.py [sys ...]"""
import subprocess, os, re, sys, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.ndimage import map_coordinates
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
QE="/workspace/qe-7.5"; CRITIC="/workspace/dimitry/calc/qe/i2_mof/tools/critic2_run.sh"
SYSTEMS=["Y_str","Tb_str","Gd_str","Y_act","Tb_act","Gd_act"]
METAL={"Y_str":"Y","Tb_str":"Tb","Gd_str":"Gd","Y_act":"Y","Y_actc":"Y","Tb_act":"Tb","Gd_act":"Gd","Gd_deep":"Gd","Tb_deep":"Tb","Y_deep":"Y"}
BOHR=0.5291772109
def sh(c):
    r=subprocess.run(c,shell=True,capture_output=True,text=True); return r.returncode,r.stdout,r.stderr
# ---------- cube I/O ----------
def read_cube(fn):
    L=open(fn).readlines(); nat=int(L[2].split()[0])
    nx,ny,nz=[int(L[i].split()[0]) for i in (3,4,5)]
    ga=np.array([float(x) for x in L[3].split()[1:4]]); gb=np.array([float(x) for x in L[4].split()[1:4]]); gc=np.array([float(x) for x in L[5].split()[1:4]])
    d=np.array([float(x) for l in L[6+nat:] for x in l.split()]).reshape((nx,ny,nz))
    return d,nx,ny,nz,np.array([ga*nx,gb*ny,gc*nz]).T,L[:6+nat]  # cellB(bohr) cols, hdr
def write_cube(fn,hdr,arr):
    with open(fn,"w") as o:
        o.writelines(hdr); f=arr.ravel(order='C')
        for i in range(0,len(f),6): o.write("".join(f"{v:13.5E}" for v in f[i:i+6])+"\n")
# ---------- atoms ----------
def read_atoms(sys):
    s=open(f"{BASE}/{sys}/scf_loaded.in").read()
    cellA=np.array([list(map(float,re.findall(r"-?\d+\.\d+",l))) for l in s.split("CELL_PARAMETERS")[1].splitlines()[1:4]])
    apos=[]
    for l in s.split("ATOMIC_POSITIONS")[1].split("K_POINTS")[0].splitlines():
        m=re.match(r"\s*([A-Z][a-z]?)\s+(-?\d.\d+)\s+(-?\d.\d+)\s+(-?\d.\d+)",l)
        if m: apos.append((m.group(1),np.array([float(m.group(2)),float(m.group(3)),float(m.group(4))])))
    return cellA,apos
# ---------- pp.x ----------
def pp(sys,prefix,plotnum,outfile):
    d=f"{BASE}/{sys}"
    typ=prefix.split('_')[-1]
    cri=f"&INPUTPP\n  prefix = '{prefix}'\n  outdir = './tmp_{typ}'\n  plot_num = {plotnum}\n  filplot = 'f'\n/\n&PLOT\n  nfile = 1\n  filepp(1) = 'f'\n  weight(1) = 1.0\n  iflag = 3\n  output_format = 6\n  fileout = '{outfile}'\n/\n"
    open(f"{d}/pp.in","w").write(cri)
    r,_,e=sh(f"cd {d} && CUDA_VISIBLE_DEVICES=0 {QE}/bin/pp.x -inp pp.in > pp_{outfile.replace('.cube','')}.log 2>&1")
    return os.path.exists(f"{d}/{outfile}")
# ---------- critic2 BCP ----------
def run_bcp(sys):
    d=f"{BASE}/{sys}"
    open(f"{d}/bcp.cri","w").write("crystal rho17.cube\nload rho17.cube id rho\nauto\n")
    sh(f"cd {d} && {CRITIC} bcp.cri > bcp.out 2>&1")
def parse_bcp(sys):
    t=open(f"{BASE}/{sys}/bcp.out",errors="ignore").read()
    cellA,apos=read_atoms(sys); cellT=cellA.T  # columns a,b,c (Å); cellT@frac=cartesian
    def near2(fc):
        res=[]
        for sym,f in apos:
            df=f-fc; df=df-np.round(df)        # fractional-space mic
            dA=cellT@df                         # cartesian displacement (Å)
            res.append((np.linalg.norm(dA),sym))
        res.sort(); return res[0][1],res[1][1],res[0][0],res[1][0]
    rows=[]
    for l in t.splitlines():
        m=re.match(r"\s*\d+\s+\w+\s+\(3,-1\)\s+\S+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+\d+\s+\w+\s+(-?\d+\.\d+E[+-]\d+)\s+\S+\s+(-?\d+\.\d+E[+-]\d+)",l)
        if m:
            fc=np.array([float(m.group(1)),float(m.group(2)),float(m.group(3))]); rho=float(m.group(4)); lap=float(m.group(5))
            a1,a2,r1,r2=near2(fc); lab="".join(sorted([a1,a2]))
            minr=min(r1,r2)
            if minr<0.55 or rho>0.15 or (r1+r2)>4.5: continue  # near-nuclear / artifact BCPs
            G=2.871*rho**(5/3)+lap/6.0  # Abramov, Ha/bohr3
            V=0.25*lap-2*G; H=G+V; De=0.5*abs(V)*2625.5  # kJ/mol
            rows.append((lab,r1+r2,rho,lap,G,V,H,De))
    return rows
# ---------- critic2 Bader CT ----------
def run_bader_ct(sys):
    d=f"{BASE}/{sys}"
    open(f"{d}/ct.cri","w").write("crystal rho17.cube\nload rho17.cube id rho\nload deltarho.cube id dr\nintegrable dr\nbader\n")
    sh(f"cd {d} && {CRITIC} ct.cri > ct.out 2>&1")
def parse_ct(sys):
    # critic2 "Integrated atomic properties" table:
    # Id cp ncp Name Z mult Volume Pop Lap $dr   ($dr = basin integral of deltarho)
    t=open(f"{BASE}/{sys}/ct.out",errors="ignore").read()
    from collections import defaultdict
    q=defaultdict(float); n=defaultdict(int)
    blk=t.split("Integrated atomic properties")[-1] if "Integrated atomic properties" in t else t
    for l in blk.splitlines():
        c=l.split()
        if len(c)>=10 and re.match(r"^\d+$",c[0]) and c[3] in ("C","H","O","Y","Tb","Gd","I"):
            try:
                dr=float(c[-1]); sym=c[3]; q[sym]+=dr; n[sym]+=1
            except ValueError: pass
    return q,n
# ---------- ELF slice ----------
def elf_slice(sys):
    d=f"{BASE}/{sys}"
    elf,nx,ny,nz,cellB,_=read_cube(f"{d}/elf.cube"); cellA=cellB.T*BOHR
    cellA=np.array([cellA[i] for i in range(3)])  # rows a,b,c in Å
    _,apos=read_atoms(sys)
    chain=np.array([f for s,f in apos if s=='I' and f[2]>0.5])  # upper-z chain group
    if len(chain)<3: chain=np.array([f for s,f in apos if s=='I'])[:4]
    cartA=(cellA@chain.T).T; cenA=cartA.mean(0)
    _,_,vh=np.linalg.svd(cartA-cenA); e1,e2,nrm=vh[0],vh[1],vh[2]; cenB=cenA/BOHR
    invB=np.linalg.inv(cellB)
    R=5.5; N=600; sg=np.linspace(-R,R,N)/BOHR; S,T=np.meshgrid(sg,sg)
    pts=cenB[None,None,:]+S[:,:,None]*e1[None,None,:]+T[:,:,None]*e2[None,None,:]
    fr=np.einsum('ij,kli->klj',invB,pts); fr-=np.floor(fr)
    vals=map_coordinates(elf,[fr[:,:,0].ravel()*nx,fr[:,:,1].ravel()*ny,fr[:,:,2].ravel()*nz],order=1,mode='wrap').reshape((N,N))
    atlist=[]
    for sym,f in apos:
        c=cellA@f; dd=c-cenA; atlist.append((sym,dd@e1,dd@e2,abs(dd@nrm)))
    fig,ax=plt.subplots(1,2,figsize=(12,5.8))
    for a,vm in zip(ax,(0.0,0.5)):
        im=a.imshow(vals,origin='lower',extent=[-R,R,-R,R],cmap='magma',vmin=vm,vmax=1.0,interpolation='bilinear')
        a.set_aspect('equal'); a.set_xlabel('e1 along chain (Å)'); a.set_ylabel('e2 (Å)')
        for sym,s,t,h in atlist:
            if abs(s)<R and abs(t)<R and h<0.5:
                a.scatter([s],[t],s=70 if sym=='I' else 25,facecolor='none',
                          edgecolor='lime' if sym=='I' else('red' if sym=='O' else'white'),lw=1.6 if sym=='I' else 1.0,zorder=5)
                if sym=='I': a.text(s,t,'  I',color='white',fontsize=8,fontweight='bold',va='center')
        plt.colorbar(im,ax=a,fraction=0.046,label='ELF')
    plt.tight_layout(); plt.savefig(f"{d}/{sys}_elf.png",dpi=150,bbox_inches='tight'); plt.close()
# ---------- driver ----------
def analyze(sys):
    d=f"{BASE}/{sys}"; print(f"\n=== {sys} ===",flush=True)
    m=METAL[sys]
    # cubes
    for pf,pn,out in [(f"{sys}_loaded",17,"rho17.cube"),(f"{sys}_loaded",8,"elf.cube"),
                      (f"{sys}_loaded",0,"val.cube"),(f"{sys}_empty",0,"empty_val.cube"),(f"{sys}_i2only",0,"i2only_val.cube")]:
        if not os.path.exists(f"{d}/{out}"):
            ok=pp(sys,pf,pn,out); print(f"  pp {out}: {'ok' if ok else 'FAIL'}",flush=True)
    # CT deltarho
    if not os.path.exists(f"{d}/deltarho.cube"):
        lv,_,_,_,_,lh=read_cube(f"{d}/val.cube"); ev,_,_,_,_,_=read_cube(f"{d}/empty_val.cube"); iv,_,_,_,_,_=read_cube(f"{d}/i2only_val.cube")
        write_cube(f"{d}/deltarho.cube",lh,lv-ev-iv)
        print("  deltarho.cube written",flush=True)
    # BCP + EML
    if not os.path.exists(f"{d}/bcp.out"): run_bcp(sys)
    bonds=parse_bcp(sys)
    # Bader CT
    if not os.path.exists(f"{d}/ct.out"): run_bader_ct(sys)
    q,n=parse_ct(sys)
    # ELF slice
    elf_slice(sys)
    # write tables
    with open(f"{d}/{sys}_bcp_eml.txt","w") as f:
        f.write(f"{sys} bond CPs + EML (Abramov G, V=¼∇²ρ−2G, D_e=½|V|)\n")
        f.write(f"{'bond':<5}{'d(Å)':>6}{'ρ':>9}{'∇²ρ':>9}{'G(Ha/b³)':>9}{'V(Ha/b³)':>9}{'H':>9}{'D_e kJ':>8}\n")
        for lab,dist,rho,lap,G,V,H,De in sorted(bonds,key=lambda x:-x[7]):
            f.write(f"{lab:<5}{dist:>6.2f}{rho:>9.4f}{lap:>9.4f}{G:>9.4f}{V:>9.4f}{H:>9.4f}{De:>8.1f}\n")
    with open(f"{d}/{sys}_ct.txt","w") as f:
        f.write(f"{sys} CT (Bader-basin δρ integral, e)\n")
        f.write(f"{'spec':<5}{'N':>4}{'Σδρ(e)':>10}{'per-atom':>10}\n")
        for s_ in ['C','H','O',m,'I']:
            if n.get(s_): f.write(f"{s_:<5}{n[s_]:>4}{q[s_]:>+10.3f}{q[s_]/n[s_]:>+10.3f}\n")
    open(f"{d}/.analyzed","w").write("ok")
    print(f"  bonds: {len(bonds)} ; wrote _bcp_eml.txt, _ct.txt, {sys}_elf.png",flush=True)
    return bonds,q,n

if __name__=="__main__":
    todo=sys.argv[1:] or SYSTEMS
    alldone=all(os.path.exists(f"{BASE}/{s}/.analyzed") for s in todo)
    for s in todo:
        if os.path.exists(f"{BASE}/{s}/.analyzed"): print(f"{s}: already analyzed, skip"); continue
        if not all(os.path.exists(f"{BASE}/{s}/{c}.done") for c in ("loaded","empty","i2only")):
            print(f"{s}: SCFs not all done, skip"); continue
        analyze(s)
