#!/usr/bin/env python3
"""Density k-mesh test on Y_deep (narrowest band). Run SCF(60/800/1e-8) then
pp.x(plot_num=17 PAW valence) for 2x2x1 (GPU0) and 4x4x2 (GPU1). Existing Γ
rho17.cube already copied as rho17_gamma.cube for comparison. Touches nothing
outside kmesh_test/."""
import subprocess,time,os
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
env=os.environ.copy()
env["PATH"]=f"{MPI}/bin:{NVHPC}/compilers/bin:{env.get('PATH','')}"
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/13.0/lib64:{env.get('LD_LIBRARY_PATH','')}"
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/kmesh_test"
def run(gpu,tag):
    log=open(f"{T}/run_{tag}.log","w")
    def sh(cmd,msg):
        log.write(msg+"\n"); log.flush(); print(f"[{time.strftime('%H:%M:%S')}] {tag}: {msg}",flush=True)
        r=subprocess.run(cmd,shell=True,env=env)
        return r.returncode
    # SCF
    rc=sh(f"cd {T} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp scf_{tag}.in > scf_{tag}.out 2>&1","SCF start")
    conv = rc==0 and "convergence has been achieved" in open(f"{T}/scf_{tag}.out",errors="ignore").read()
    sh("(scf done)" if conv else f"(SCF rc={rc}!)","SCF "+("CONVERGED" if conv else "FAILED"))
    if not conv: log.close(); return
    # pp plot_num=17
    rc=sh(f"cd {T} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pp.x -inp pp_{tag}.in > pp_{tag}.out 2>&1","pp(17) start")
    ok = os.path.exists(f"{T}/rho17_{tag}.cube") and os.path.getsize(f"{T}/rho17_{tag}.cube")>1e7
    sh("(pp done)" if ok else f"(pp rc={rc}!)","pp(17) "+("CUBE WRITTEN" if ok else "FAILED"))
    log.close()
import threading
t1=threading.Thread(target=run,args=(0,"k221"))
t2=threading.Thread(target=run,args=(1,"k442"))
t1.start(); t2.start(); t1.join(); t2.join()
print("=== kmesh density test COMPLETE ===",flush=True)
