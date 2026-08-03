#!/usr/bin/env python3
"""Run the 4 finite-displacement SCFs (gp/gm Gamma +/-, k442p/k442m 4x4x2 +/-) for
the I2-I3 stretch force-constant spot-check. 2 GPUs, memory-checked."""
import subprocess,time,os,glob
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
env=os.environ.copy()
env["PATH"]=f"{MPI}/bin:{NVHPC}/compilers/bin:{env.get('PATH','')}"
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/13.0/lib64:{env.get('LD_LIBRARY_PATH','')}"
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/kmesh_test/raman_fd"
def gpuram(i): return int(subprocess.run(f"nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i {i}",shell=True,capture_output=True,text=True).stdout.strip() or 99999)
JOBS=["gp","gm","k442p","k442m"]; q=list(JOBS); running={}
def launch(gpu,tag):
    for tmp in glob.glob(f"{T}/tmp_{tag}*"): subprocess.run(f"rm -rf {tmp}",shell=True)
    cmd=f"cd {T} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp scf_{tag}.in > scf_{tag}.out 2>&1"
    return subprocess.Popen(cmd,shell=True,env=env)
def trylaunch(gpu):
    if gpu in running or not q: return
    if gpuram(gpu)>2000: return
    tag=q.pop(0); running[gpu]=(tag,launch(gpu,tag),time.time())
    print(f"[{time.strftime('%H:%M:%S')}] GPU{gpu} START {tag}",flush=True)
print("4 FD SCFs (Gamma fast, 4x4x2 slower)",flush=True)
while running or q:
    for g in (0,1): trylaunch(g)
    time.sleep(15)
    for g in list(running):
        tag,p,t0=running[g]
        if p.poll() is not None:
            ok="convergence has been achieved" in open(f"{T}/scf_{tag}.out",errors="ignore").read()
            print(f"[{time.strftime('%H:%M:%S')}] GPU{g} DONE {tag} rc={p.returncode} conv={ok} ({time.time()-t0:.0f}s)",flush=True)
            del running[g]
print("=== 4 FD SCFs COMPLETE ===",flush=True)
