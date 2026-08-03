#!/usr/bin/env python3
import subprocess,time,os,glob
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
env=os.environ.copy()
env["PATH"]=f"{MPI}/bin:{NVHPC}/compilers/bin:{env.get('PATH','')}"
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/13.0/lib64:{env.get('LD_LIBRARY_PATH','')}"
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/relax_442/sp02"
QUEUE=["Gdact","Tbact","Yst","Gdst","Tbst"]
def gpuram(i): return int(subprocess.run(f"nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i {i}",shell=True,capture_output=True,text=True).stdout.strip() or 99999)
def done(m): return os.path.exists(f"{T}/sp02_{m}.out") and "convergence has been achieved" in open(f"{T}/sp02_{m}.out",errors="ignore").read() and "!" in open(f"{T}/sp02_{m}.out",errors="ignore").read()
q=[m for m in QUEUE if not done(m)]; running={}; log=open(f"{T}/run.log","w")
def launch(g,m): 
    subprocess.run(f"rm -rf {T}/tmp_sp02_{m}*",shell=True)
    cmd=f"cd {T} && CUDA_VISIBLE_DEVICES={g} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp sp02_{m}.in > sp02_{m}.out 2>&1"
    return subprocess.Popen(cmd,shell=True,env=env,start_new_session=True)
while running or q:
    for g in (0,1):
        if g in running or not q: continue
        if gpuram(g)>3000: continue
        m=q.pop(0); running[g]=(m,launch(g,m)); log.write(f"[{time.strftime('%H:%M:%S')}] GPU{g} SP {m}\n"); log.flush()
    time.sleep(30)
    for g in list(running):
        m,p=running[g]
        if p.poll() is not None: log.write(f"[{time.strftime('%H:%M:%S')}] GPU{g} DONE {m} ok={done(m)}\n"); log.flush(); del running[g]
log.write("=== SP02 complete ===\n")
