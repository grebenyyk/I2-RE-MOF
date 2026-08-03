#!/usr/bin/env python3
"""Driver: run 9 iodine-only 4×4×2 relaxations, 2 at a time (1/GPU).
Memory-checked. Sorted fast→slow (vc first — least metallic, easiest SCF)."""
import subprocess,time,os,glob
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
env=os.environ.copy()
env["PATH"]=f"{MPI}/bin:{NVHPC}/compilers/bin:{env.get('PATH','')}"
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/13.0/lib64:{env.get('LD_LIBRARY_PATH','')}"
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/relax_442"
# order: vc first (least metallic), then act, then st
ORDER=["Yvc","Gdvc","Tbvc","Yact","Gdact","Tbact","Yst","Gdst","Tbst"]
JOBS=[f"relax442_{s}" for s in ORDER]
def gpuram(i): return int(subprocess.run(f"nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i {i}",shell=True,capture_output=True,text=True).stdout.strip() or 99999)
def converged(tag):
    o=f"{T}/{tag}.out"
    return os.path.exists(o) and "bfgs converged" in open(o,errors="ignore").read()
todo=[t for t in JOBS if not converged(t)]
log=open(f"{T}/run.log","w"); running={}; q=list(todo); done=0; tot=len(todo)
def launch(gpu,tag):
    for tmp in glob.glob(f"{T}/tmp_{tag}*"): subprocess.run(f"rm -rf {tmp}",shell=True)
    cmd=f"cd {T} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp {tag}.in > {tag}.out 2>&1"
    return subprocess.Popen(cmd,shell=True,env=env)
def trylaunch(gpu):
    if gpu in running or not q: return
    if gpuram(gpu)>3000: return
    tag=q.pop(0); running[gpu]=(tag,launch(gpu,tag),time.time())
    m=f"[{time.strftime('%H:%M:%S')}] GPU{gpu} START {tag} [{tot-len(q)}/{tot}]"; print(m,flush=True); log.write(m+"\n"); log.flush()
print(f"{tot} 4×4×2 iodine relaxations to run (vc→act→st)",flush=True)
while running or q:
    for g in (0,1): trylaunch(g)
    time.sleep(20)
    for g in list(running):
        tag,p,t0=running[g]
        if p.poll() is not None:
            ok=converged(tag); done+=1
            e="";
            try:
                import re; m=re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f"{T}/{tag}.out",errors="ignore").read(),re.M); e=m[-1] if m else "?"
            except: pass
            m=f"[{time.strftime('%H:%M:%S')}] GPU{g} DONE {tag} conv={ok} E={e} ({time.time()-t0:.0f}s) [{done}/{tot}]"; print(m,flush=True); log.write(m+"\n"); log.flush()
            del running[g]
print(f"=== ALL {done} RELAXATIONS COMPLETE ===",flush=True)
