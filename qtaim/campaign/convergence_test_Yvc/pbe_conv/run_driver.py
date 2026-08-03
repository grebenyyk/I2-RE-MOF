#!/usr/bin/env python3
"""Driver: run the 9 PBE+D3 convergence SCFs 2 at a time (1 per L40S GPU).
Queues jobs; when a GPU frees up, launches the next. Restarts are clean
(old .out removed). Matches the campaign convention (NP=1, 1 GPU per calc)."""
import subprocess, time, os, glob
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/pbe_conv"
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"
QE="/workspace/qe-7.5"
env=os.environ.copy()
env["PATH"]=MPI["PATH"] if False else f"{MPI}/bin:{NVHPC}/compilers/bin:{env.get('PATH','')}"
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/13.0/lib64:{env.get('LD_LIBRARY_PATH','')}"
JOBS=["pbe_gas_ec40","pbe_gas_ec60","pbe_empty_ec40","pbe_empty_ec60",
      "pbe_kG","pbe_k221","pbe_poly_ec40","pbe_poly_ec60","pbe_k663"]
running={}  # gpu -> (tag, popen, t0)
done=set()
log=open(f"{D}/driver.log","w")
def launch(gpu):
    if not JOBS: return
    tag=JOBS.pop(0)
    # clean old outputs
    for p in glob.glob(f"{D}/tmp_{tag}*"): subprocess.run(f"rm -rf {p}",shell=True)
    cmd=f"cd {D} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp {tag}.in > {tag}.out 2>&1"
    running[gpu]=(tag,subprocess.Popen(cmd,shell=True,env=env),time.time())
    msg=f"[{time.strftime('%H:%M:%S')}] GPU{gpu} START {tag}"
    print(msg,flush=True); log.write(msg+"\n"); log.flush()
print(f"{len(JOBS)} PBE+D3 convergence calcs to run (2 GPUs, 1 calc/GPU)",flush=True)
# launch first 2
launch(0); launch(1)
while running:
    time.sleep(20)
    for gpu in list(running):
        tag,p,t0=running[gpu]
        if p.poll() is not None:
            dt=time.time()-t0
            ok="convergence has been achieved" in open(f"{D}/{tag}.out",errors="ignore").read()
            msg=f"[{time.strftime('%H:%M:%S')}] GPU{gpu} DONE {tag} rc={p.returncode} conv={ok} ({dt:.0f}s)"
            print(msg,flush=True); log.write(msg+"\n"); log.flush()
            done.add(tag); del running[gpu]
            if JOBS: launch(gpu)
print(f"=== ALL {len(done)} CONVERGENCE CALCS COMPLETE ===",flush=True)
for tag in ["pbe_kG","pbe_k221","pbe_k663","pbe_poly_ec40","pbe_poly_ec60",
            "pbe_empty_ec40","pbe_empty_ec60","pbe_gas_ec40","pbe_gas_ec60"]:
    e=grep=f"  {tag}: " + subprocess.run(f"grep -aE '!.*total energy' {D}/{tag}.out 2>/dev/null | tail -1",
                                         shell=True,capture_output=True,text=True).stdout.strip()
    print(e,flush=True)
