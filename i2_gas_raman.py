#!/usr/bin/env python3
"""Run gas-phase I2 through the SAME finite-difference pipeline as I2@Y-MOF, to
validate the pipeline and the PBE/PAW level. Expected: I-I stretch ~213 cm-1.
If it comes out high (e.g. >225), the pipeline/level carries a systematic
blueshift and the MOF Raman should be presented as fingerprints, not absolutes.

Mirrors raman_driver.py: delta=0.01 A, central differences, Hessian in Ry/bohr^2,
mass-weight (M_I=126.904), RY_TO_CM=5140.4843. 2 I atoms in a 15 A vacuum box;
relax first, then 12 displaced SCFs, 6x6 Hessian -> frequencies.
ecutrho 400 here (vs 800 in MOF): frequencies are insensitive (<~2 cm-1)."""
import os, re, subprocess, numpy as np
env=os.environ.copy()
MPI="/workspace/openmpi"; NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; CUDA="13.0"
env["PATH"]=f"{MPI}/bin:"+env["PATH"]
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/{CUDA}/lib64:{NVHPC}/cuda/{CUDA}/lib64"
env["OMP_NUM_THREADS"]="8"
QE="/workspace/qe-7.5/bin/pw.x"
D="/workspace/dimitry/calc/qe/i2_mof/i2_gas_raman"; os.makedirs(D, exist_ok=True)
ECUTWFC=60; ECUTRHO=400; ABOX=15.0; R0=2.67
DELTA_A=0.01; M_I=126.904; RY_TO_CM=5140.4843; BOHR=0.529177210563

def write_scf(path,prefix,pos,calc="scf",outdir="./tmp"):
    L=[f"&CONTROL\n  calculation = '{calc}'\n  prefix = '{prefix}'\n  outdir = '{outdir}'\n",
       "  pseudo_dir = '/workspace/qe-7.5/pseudo'\n  tprnfor = .true.\n/\n",
       "&SYSTEM\n  ibrav = 0\n  nat = 2\n  ntyp = 1\n",
       f"  ecutwfc = {ECUTWFC}\n  ecutrho = {ECUTRHO}\n  occupations = 'fixed'\n  nspin = 1\n/\n",
       "&ELECTRONS\n  conv_thr = 1.0d-8\n  diagonalization = 'david'\n  mixing_beta = 0.2\n/\n",
       "&IONS\n/\n",
       "ATOMIC_SPECIES\n  I  126.904  I.pbe-n-kjpaw_psl.1.0.0.UPF\n",
       "CELL_PARAMETERS angstrom\n",
       f"  {ABOX} 0 0\n  0 {ABOX} 0\n  0 0 {ABOX}\n","ATOMIC_POSITIONS angstrom\n"]
    for p in pos: L.append(f"  I  {p[0]:.10f}  {p[1]:.10f}  {p[2]:.10f}\n")
    L.append("K_POINTS gamma\n"); open(path,"w").write("".join(L))

def run(inp,out,outdir):
    subprocess.run(f"cd {D} && rm -rf {outdir} && CUDA_VISIBLE_DEVICES=0 mpirun --allow-run-as-root -np 1 {QE} -nk 1 -nd 1 -nb 1 -nt 1 -inp {inp} > {out} 2>&1",
                   shell=True, env=env, check=True)

def parse_atoms_ang(t):
    blk=re.split(r"ATOMIC_POSITIONS \(angstrom\)",t)[-1].split("K_POINTS")[0]
    return np.array([[float(x) for x in re.findall(r"-?\d+\.\d+",ln)]
                     for ln in blk.split("\n") if re.match(r"\s*I\s+-?\d",ln)])
def parse_forces(t):
    blk=re.split(r"Forces acting on atoms",t)[-1]
    f=[]
    for ln in blk.split("\n"):
        m=re.match(r"\s*atom\s+\d+\s+type\s+\S+\s+force\s+=\s+(-?\S+)\s+(-?\S+)\s+(-?\S+)",ln)
        if m: f.append([float(m.group(k)) for k in (1,2,3)])
    return np.array(f)  # (2,3) Ry/bohr

# 1. relax I2
pos0=np.array([[ABOX/2-R0/2,ABOX/2,ABOX/2],[ABOX/2+R0/2,ABOX/2,ABOX/2]])
write_scf(f"{D}/relax.in","i2gas",pos0,calc="relax",outdir="./tmp_rel")
print("relaxing I2 ..."); run("relax.in","relax.out","tmp_rel")
posr=parse_atoms_ang(open(f"{D}/relax.out").read())
bond=np.linalg.norm(posr[1]-posr[0])
print(f"relaxed I-I = {bond:.4f} A  (exp 2.67)")

# 2. finite differences (12 SCFs)
delta_bohr=DELTA_A/BOHR; F={}
for i in range(2):
    for a in range(3):
        for s,sgn in [('p',1),('n',-1)]:
            p=posr.copy(); p[i,a]+=sgn*DELTA_A
            tag=f"I{i}_{['x','y','z'][a]}{s}"
            write_scf(f"{D}/disp_{tag}.in",f"i2gas_{tag}",p,calc="scf",outdir=f"./tmp_{tag}")
            run(f"disp_{tag}.in",f"disp_{tag}.out",f"tmp_{tag}")
            F[(i,a,s)]=parse_forces(open(f"{D}/disp_{tag}.out").read())
            print(f"  {tag} done")

# 3. Hessian 6x6
Phi=np.zeros((6,6))
for i in range(2):
    for a in range(3):
        col=3*i+a
        dF=F[(i,a,'p')]-F[(i,a,'n')]
        Phi[:,col]=-dF.flatten()/(2*delta_bohr)
masses=np.array([M_I]*6)
Dm=Phi/np.sqrt(masses[:,None]*masses[None,:])
ev,vec=np.linalg.eigh(Dm)
freqs=np.sign(ev)*np.sqrt(np.abs(ev))*RY_TO_CM
order=np.argsort(freqs)[::-1]
print("\n=== gas I2 frequencies (cm-1) ===")
for k,ix in enumerate(order):
    print(f"  mode {k+1}: {freqs[ix]:8.1f}   (eig {ev[ix]:.4e})")
print(f"\nstretch (highest non-zero) ~ {freqs[order[0]]:.1f} cm-1  (exp gas I2 = 213)")
np.save(f"{D}/freqs_gas.npy",freqs); np.save(f"{D}/Phi_gas.npy",Phi)
print("DONE")
