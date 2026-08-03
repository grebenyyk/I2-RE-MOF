#!/usr/bin/env python3
"""Run empty-MOF k-mesh series (2×2×1, 4×4×2, 6×6×3) + gas 2×2×1 confirmation.
Proves whether empty MOF & gas are converged at Γ (current production assumption).
Γ empty = rd_empty_Yvc (already done, genuine PBE). 2 GPUs, memory-checked."""
import subprocess,time,os,glob
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/pbe_conv"
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
env=os.environ.copy()
env["PATH"]=f"{MPI}/bin:{NVHPC}/compilers/bin:{env.get('PATH','')}"
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/13.0/lib64:{env.get('LD_LIBRARY_PATH','')}"
JOBS=["pbe_gas_k221","pbe_empty_k221","pbe_empty_k442","pbe_empty_k663"]
def gpuram(i): return int(subprocess.run(f"nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i {i}",shell=True,capture_output=True,text=True).stdout.strip() or 99999)
log=open(f"{D}/empty_kseries.log","w"); running={}; q=list(JOBS)
def launch(gpu,tag):
    for tmp in glob.glob(f"{D}/tmp_{tag}*"): subprocess.run(f"rm -rf {tmp}",shell=True)
    cmd=f"cd {D} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp {tag}.in > {tag}.out 2>&1"
    return subprocess.Popen(cmd,shell=True,env=env)
def trylaunch(gpu):
    if gpu in running or not q: return
    if gpuram(gpu)>2000: return
    tag=q.pop(0); running[gpu]=(tag,launch(gpu,tag),time.time())
    msg=f"[{time.strftime('%H:%M:%S')}] GPU{gpu} START {tag}"; print(msg,flush=True); log.write(msg+"\n"); log.flush()
print(f"{len(JOBS)} calcs: empty k-series + gas confirmation",flush=True)
while running or q:
    for g in (0,1): trylaunch(g)
    time.sleep(15)
    for g in list(running):
        tag,p,t0=running[g]
        if p.poll() is not None:
            ok="convergence has been achieved" in open(f"{D}/{tag}.out",errors="ignore").read()
            msg=f"[{time.strftime('%H:%M:%S')}] GPU{g} DONE {tag} rc={p.returncode} conv={ok} ({time.time()-t0:.0f}s)"; print(msg,flush=True); log.write(msg+"\n"); log.flush()
            del running[g]
print("=== empty k-series + gas confirmation COMPLETE ===",flush=True)
