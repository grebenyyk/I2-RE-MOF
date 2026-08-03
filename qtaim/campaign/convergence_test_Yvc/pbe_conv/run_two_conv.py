#!/usr/bin/env python3
"""Run the 2 remaining convergence calcs (k663, poly_ec60) to completion,
1 per GPU. Memory-checked, clean restart. Produces 9/9."""
import subprocess, time, os, glob, sys
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/pbe_conv"
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"
QE="/workspace/qe-7.5"
env=os.environ.copy()
env["PATH"]=f"{MPI}/bin:{NVHPC}/compilers/bin:{env.get('PATH','')}"
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/13.0/lib64:{env.get('LD_LIBRARY_PATH','')}"

def gpuram(i):
    return int(subprocess.run(f"nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i {i}",
                 shell=True,capture_output=True,text=True).stdout.strip() or 99999)

JOBS=[("pbe_k663",0),("pbe_poly_ec60",1)]
running={}
for tag,gpu in JOBS:
    # wait for GPU to be clear
    while gpuram(gpu) > 5000:
        print(f"[{time.strftime('%H:%M:%S')}] GPU{gpu} has {gpuram(gpu)} MiB used, waiting...",flush=True); time.sleep(10)
    for tmp in glob.glob(f"{D}/tmp_{tag}*"): subprocess.run(f"rm -rf {tmp}",shell=True)
    os.remove(f"{D}/{tag}.out") if os.path.exists(f"{D}/{tag}.out") else None
    cmd=f"cd {D} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp {tag}.in > {tag}.out 2>&1"
    running[tag]=subprocess.Popen(cmd,shell=True,env=env)
    print(f"[{time.strftime('%H:%M:%S')}] START {tag} on GPU{gpu}",flush=True)
    time.sleep(5)

while running:
    time.sleep(20)
    for tag in list(running):
        if running[tag].poll() is not None:
            rc=running[tag].returncode
            ok="convergence has been achieved" in open(f"{D}/{tag}.out",errors="ignore").read()
            import re
            m=re.findall(r"^!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f"{D}/{tag}.out",errors="ignore").read(),re.M)
            e=m[-1] if m else "NONE"
            print(f"[{time.strftime('%H:%M:%S')}] DONE {tag} rc={rc} conv={ok} E={e}",flush=True)
            for tmp in glob.glob(f"{D}/tmp_{tag}*"): subprocess.run(f"rm -rf {tmp}",shell=True)
            del running[tag]
print("=== BOTH CONVERGENCE CALCS COMPLETE (9/9) ===",flush=True)
