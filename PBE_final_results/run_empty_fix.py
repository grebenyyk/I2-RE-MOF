#!/usr/bin/env python3
"""Re-run ALL empty-MOF SCFs at conv_thr=1e-6 (were 1e-3, the only spurious
inconsistency vs loaded). 20 jobs: 5 convergence + 15 production. Slow
convergence k-mesh points first. 2 GPUs, memory-checked."""
import subprocess,time,os,glob
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
env=os.environ.copy()
env["PATH"]=f"{MPI}/bin:{NVHPC}/compilers/bin:{env.get('PATH','')}"
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/13.0/lib64:{env.get('LD_LIBRARY_PATH','')}"
P="/workspace/dimitry/calc/qe/i2_mof"
PC=f"{P}/qtaim/campaign/convergence_test_Yvc/pbe_conv"
PB=f"{P}/PBE_final_results"
# (directory, tag) — slow convergence k-mesh first, then convergence Γ/ecut, then production Γ
JOBS=[
 (PC,"pbe_empty_k663"),(PC,"pbe_empty_k442"),(PC,"pbe_empty_k221"),
 (PC,"pbe_empty_ec40"),(PC,"pbe_empty_ec60"),
 (f"{PB}/adsorption_9systems","rd_empty_Yvc"),(f"{PB}/adsorption_9systems","rd_empty_Yact"),(f"{PB}/adsorption_9systems","rd_empty_Yst"),
 (f"{PB}/adsorption_9systems","rd_empty_Gdvc"),(f"{PB}/adsorption_9systems","rd_empty_Gdact"),(f"{PB}/adsorption_9systems","rd_empty_Gdst"),
 (f"{PB}/adsorption_9systems","rd_empty_Tbvc"),(f"{PB}/adsorption_9systems","rd_empty_Tbact"),(f"{PB}/adsorption_9systems","rd_empty_Tbst"),
 (f"{PB}/decomposition_Y3","rd_empty_Yvc"),(f"{PB}/decomposition_Y3","rd_empty_Yact"),(f"{PB}/decomposition_Y3","rd_empty_Yst"),
 (f"{PB}/landscape_3variants","rd_empty_Yvc"),(f"{PB}/landscape_3variants","rd_empty_Yact"),(f"{PB}/landscape_3variants","rd_empty_Yst"),
]
def gpuram(i): return int(subprocess.run(f"nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i {i}",shell=True,capture_output=True,text=True).stdout.strip() or 99999)
log=open(f"{PB}/run_empty_fix.log","w"); running={}; q=list(JOBS); done=0; total=len(JOBS)
def launch(gpu,d,tag):
    for tmp in glob.glob(f"{d}/tmp_{tag}*"): subprocess.run(f"rm -rf {tmp}",shell=True)
    cmd=f"cd {d} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp {tag}.in > {tag}.out 2>&1"
    return subprocess.Popen(cmd,shell=True,env=env)
def trylaunch(gpu):
    if gpu in running or not q: return
    if gpuram(gpu)>2000: return
    d,tag=q.pop(0); running[gpu]=(d,tag,launch(gpu,d,tag),time.time())
    msg=f"[{time.strftime('%H:%M:%S')}] GPU{gpu} START {tag} ({d.split('/PBE_final_results/')[-1] if 'PBE_final_results' in d else 'pbe_conv'})"; print(msg,flush=True); log.write(msg+"\n"); log.flush()
print(f"{total} empty re-runs at conv_thr=1e-6 (slow convergence k-mesh first)",flush=True)
while running or q:
    for g in (0,1): trylaunch(g)
    time.sleep(15)
    for g in list(running):
        d,tag,p,t0=running[g]
        if p.poll() is not None:
            ok="convergence has been achieved" in open(f"{d}/{tag}.out",errors="ignore").read()
            done+=1
            msg=f"[{time.strftime('%H:%M:%S')}] GPU{g} DONE {tag} rc={p.returncode} conv={ok} ({time.time()-t0:.0f}s) [{done}/{total}]"; print(msg,flush=True); log.write(msg+"\n"); log.flush()
            del running[g]
print(f"=== ALL {done} empty re-runs (1e-6) COMPLETE ===",flush=True)
