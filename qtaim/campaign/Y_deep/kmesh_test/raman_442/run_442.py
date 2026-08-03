#!/usr/bin/env python3
"""Run the 48 displaced 4x4x2 SCFs for the I-sublattice Hessian (2 GPUs,
memory-checked). Each seeded with Gamma charge density. ~6h."""
import subprocess,time,os,glob
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
env=os.environ.copy()
env["PATH"]=f"{MPI}/bin:{NVHPC}/compilers/bin:{env.get('PATH','')}"
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/13.0/lib64:{env.get('LD_LIBRARY_PATH','')}"
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/kmesh_test/raman_442"
tags=[l.split()[3] for l in open(f"{T}/tags.txt") if l.strip()]
def gpuram(i): return int(subprocess.run(f"nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i {i}",shell=True,capture_output=True,text=True).stdout.strip() or 99999)
def done(tag):
    o=f"{T}/disp_{tag}.out"
    return os.path.exists(o) and "Forces acting on atoms" in open(o,errors="ignore").read() and "convergence has been achieved" in open(o,errors="ignore").read()
todo=[t for t in tags if not done(t)]
log=open(f"{T}/run.log","w"); running={}; q=list(todo); n=0; tot=len(todo)
def launch(gpu,tag):
    cmd=f"cd {T} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp disp_{tag}.in > disp_{tag}.out 2>&1"
    return subprocess.Popen(cmd,shell=True,env=env)
def trylaunch(gpu):
    if gpu in running or not q: return
    if gpuram(gpu)>3000: return
    tag=q.pop(0); running[gpu]=(tag,launch(gpu,tag),time.time())
    m=f"[{time.strftime('%H:%M:%S')}] GPU{gpu} START {tag} [{tot-len(q)}/{tot}]"; print(m,flush=True); log.write(m+"\n"); log.flush()
print(f"{tot} displaced 4x4x2 SCFs to run (2 GPUs)",flush=True)
while running or q:
    for g in (0,1): trylaunch(g)
    time.sleep(20)
    for g in list(running):
        tag,p,t0=running[g]
        if p.poll() is not None:
            ok=done(tag); n+=1
            m=f"[{time.strftime('%H:%M:%S')}] GPU{g} DONE {tag} {'OK' if ok else 'CHECK'} ({time.time()-t0:.0f}s) [{n}/{tot}]"; print(m,flush=True); log.write(m+"\n"); log.flush()
            del running[g]
print(f"=== all {n} done ===",flush=True)
