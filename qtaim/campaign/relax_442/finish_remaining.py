#!/usr/bin/env python3
"""Targeted driver for the 3 remaining 4x4x2 iodine relaxations:
   Gdst, Tbst (standard seeds, bfgs) and Gdact_yseed (Y-act I seed, bfgs).
   Memory-guarded: will NOT launch on a GPU with >3 GB used, so it cannot
   disturb Yst (GPU0) or Tbact_yseed (GPU1) while they run. Picks up each
   pending calc as a GPU frees. Survives relaunch (idempotent on convergence)."""
import subprocess,time,os,glob
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
env=os.environ.copy()
env["PATH"]=f"{MPI}/bin:{NVHPC}/compilers/bin:{env.get('PATH','')}"
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/13.0/lib64:{env.get('LD_LIBRARY_PATH','')}"
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/relax_442"
QUEUE=["relax442_Gdst","relax442_Tbst","relax442_Gdact_yseed"]
def gpuram(i): return int(subprocess.run(f"nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i {i}",shell=True,capture_output=True,text=True).stdout.strip() or 99999)
def converged(tag):
    o=f"{T}/{tag}.out"
    return os.path.exists(o) and "bfgs converged" in open(o,errors="ignore").read()
def running_pwxs():
    out=subprocess.run("pgrep -af 'pw.x.*-inp relax442'",shell=True,capture_output=True,text=True).stdout
    return out
def launch(gpu,tag):
    for tmp in glob.glob(f"{T}/tmp_{tag}*"): subprocess.run(f"rm -rf {tmp}",shell=True)
    cmd=f"cd {T} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp {tag}.in > {tag}.out 2>&1"
    return subprocess.Popen(cmd,shell=True,env=env,start_new_session=True)
q=[t for t in QUEUE if not converged(t)]
log=open(f"{T}/finish.log","a"); running={}; done=0; tot=len(q)
def logline(m): print(m,flush=True); log.write(m+"\n"); log.flush()
logline(f"[{time.strftime('%H:%M:%S')}] finish_remaining start: {tot} pending {q}")
while running or q:
    for g in (0,1):
        if g in running or not q: continue
        if gpuram(g)>3000: continue          # GPU busy (Yst/Tbact_yseed or other)
        tag=q.pop(0); running[g]=(tag,launch(g,tag),time.time())
        logline(f"[{time.strftime('%H:%M:%S')}] GPU{g} START {tag} [{tot-len(q)}/{tot}]")
    time.sleep(30)
    for g in list(running):
        tag,p,t0=running[g]
        if p.poll() is not None:
            ok=converged(tag); done+=1
            try:
                import re; m=re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f"{T}/{tag}.out",errors="ignore").read(),re.M); e=m[-1] if m else "?"
            except: e="?"
            logline(f"[{time.strftime('%H:%M:%S')}] GPU{g} DONE {tag} conv={ok} E={e} ({time.time()-t0:.0f}s) [{done}/{tot}]")
            del running[g]
logline(f"[{time.strftime('%H:%M:%S')}] === finish_remaining complete: {done} done ===")
