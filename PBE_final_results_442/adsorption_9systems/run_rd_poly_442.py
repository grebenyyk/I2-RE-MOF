#!/usr/bin/env python3
"""Run the real 0.02 scf single-points (rd_poly_442) for vc + Yact on their 4x4x2-relaxed
geometries, so rd_poly_442_<sys>.in/.out are REAL corresponding files (not synthetic).
vc/Yact were relaxed at 0.02 already; this re-evaluates the energy as a clean dedicated scf.
Memory-guarded (2 GPUs), queues behind other calcs. Idempotent (skips converged)."""
import subprocess,time,os,re,fcntl,sys
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
env=os.environ.copy()
env["PATH"]=f"{MPI}/bin:{NVHPC}/compilers/bin:{env.get('PATH','')}"
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/13.0/lib64:{env.get('LD_LIBRARY_PATH','')}"
D="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results_442/adsorption_9systems"
TAGS=["Yvc","Gdvc","Tbvc","Yact"]
lf=open(f"{D}/rdpoly.lock","w")
try: fcntl.flock(lf,fcntl.LOCK_EX|fcntl.LOCK_NB)
except (BlockingIOError,OSError): sys.exit(0)
def gpuram(i): return int(subprocess.run(f"nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i {i}",shell=True,capture_output=True,text=True).stdout.strip() or 99999)
def done(t):
    o=f"{D}/rd_poly_442_{t}.out"
    return os.path.exists(o) and "convergence has been achieved" in open(o,errors="ignore").read() and "!" in open(o,errors="ignore").read()
def launch(g,t):
    subprocess.run(f"rm -rf {D}/tmp_rd_poly_442_{t}*",shell=True)
    cmd=f"cd {D} && CUDA_VISIBLE_DEVICES={g} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp rd_poly_442_{t}.in > rd_poly_442_{t}.out 2>&1"
    return subprocess.Popen(cmd,shell=True,env=env,start_new_session=True)
log=open(f"{D}/rdpoly_run.log","a")
def logl(m): log.write(f"[{time.strftime('%H:%M:%S')}] {m}\n"); log.flush()
q=[t for t in TAGS if not done(t)]; running={}
while running or q:
    for g in (0,1):
        if g in running or not q or gpuram(g)>3000: continue
        t=q.pop(0); running[g]=(t,launch(g,t)); logl(f"GPU{g} START rd_poly_442_{t}")
    time.sleep(20)
    for g in list(running):
        t,p=running[g]
        if p.poll() is not None: logl(f"GPU{g} DONE {t} ok={done(t)}"); del running[g]
logl("=== rd_poly_442 vc/Yact complete ===")
