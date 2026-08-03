#!/usr/bin/env python3
"""Y_actc density k-mesh test (wider band -> expect smaller Gamma discrepancy).
SCF(60/800/1e-8) then pp(17) for k221 and k442, memory-checked (uses free GPU).
Touches nothing outside Y_actc/kmesh_test/."""
import subprocess,time,os,glob
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
env=os.environ.copy()
env["PATH"]=f"{MPI}/bin:{NVHPC}/compilers/bin:{env.get('PATH','')}"
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/13.0/lib64:{env.get('LD_LIBRARY_PATH','')}"
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_actc/kmesh_test"
def gpuram(i): return int(subprocess.run(f"nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i {i}",shell=True,capture_output=True,text=True).stdout.strip() or 99999)
log=open(f"{T}/run.log","w")
JOBS=["k221","k442"]; q=list(JOBS); running={}
def run_scf_pp(gpu,tag):
    for tmp in glob.glob(f"{T}/tmp_{tag}*"): subprocess.run(f"rm -rf {tmp}",shell=True)
    sc=f"cd {T} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp scf_{tag}.in > scf_{tag}.out 2>&1"
    r=subprocess.run(sc,shell=True,env=env)
    conv = r.returncode==0 and "convergence has been achieved" in open(f"{T}/scf_{tag}.out",errors="ignore").read()
    msg=f"[{time.strftime('%H:%M:%S')}] GPU{gpu} {tag} SCF {'CONVERGED' if conv else 'FAILED rc='+str(r.returncode)}"; print(msg,flush=True); log.write(msg+"\n"); log.flush()
    if not conv: return
    pp=f"cd {T} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pp.x -inp pp_{tag}.in > pp_{tag}.out 2>&1"
    subprocess.run(pp,shell=True,env=env)
    ok=os.path.exists(f"{T}/rho17_{tag}.cube") and os.path.getsize(f"{T}/rho17_{tag}.cube")>1e7
    msg=f"[{time.strftime('%H:%M:%S')}] GPU{gpu} {tag} pp(17) {'CUBE WRITTEN' if ok else 'FAILED'}"; print(msg,flush=True); log.write(msg+"\n"); log.flush()
import threading
def worker(gpu,tag):
    # wait for GPU free
    while gpuram(gpu)>2000: time.sleep(15)
    run_scf_pp(gpu,tag)
# assign: k221 -> GPU1 (free now), k442 -> GPU0 (when k663 done)
t1=threading.Thread(target=worker,args=(1,"k221")); t2=threading.Thread(target=worker,args=(0,"k442"))
t1.start(); t2.start(); t1.join(); t2.join()
print("=== Y_actc kmesh density test COMPLETE ===",flush=True)
