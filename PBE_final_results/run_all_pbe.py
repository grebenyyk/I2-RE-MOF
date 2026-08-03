#!/usr/bin/env python3
"""Driver: regenerate ALL genuine-PBE+D3 SCFs in PBE_final_results/.
- PRIORITY: the 5 Y-decomposition pieces first (so the Y table is ready ASAP):
    decomposition_Y3/rd_poly_Y{st,act,vc}, rd_empty_Y{act,vc}
- 1 calc per GPU; launches on a GPU only when it is essentially free (<2 GB used),
  so it coexists with a still-running convergence calc without OOM.
- Skips any .in whose .out is already converged AND genuine PBE
  (discriminator "Exchange-correlation= PBE"; HSE .out says "= HSE" → not skipped).
"""
import subprocess, time, os, glob
PBE="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results"
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"
QE="/workspace/qe-7.5"
env=os.environ.copy()
env["PATH"]=f"{MPI}/bin:{NVHPC}/compilers/bin:{env.get('PATH','')}"
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/13.0/lib64:{env.get('LD_LIBRARY_PATH','')}"

def is_pbe_done(inp):
    out=inp[:-3]+".out"
    if not os.path.exists(out): return False
    txt=open(out,errors="ignore").read()
    return ("convergence has been achieved" in txt) and ("Exchange-correlation= PBE" in txt)

all_ins=[f for f in glob.glob(f"{PBE}/*/*.in") if not is_pbe_done(f)]
PRIORITY={"/decomposition_Y3/rd_poly_Yst.in","/decomposition_Y3/rd_poly_Yact.in",
          "/decomposition_Y3/rd_poly_Yvc.in","/decomposition_Y3/rd_empty_Yact.in",
          "/decomposition_Y3/rd_empty_Yvc.in"}
def prio(f):
    return (0 if f.endswith(tuple(PRIORITY)) else (1 if "/decomposition_Y3/" in f else 2),
            os.path.getsize(f))
all_ins.sort(key=prio)

def gpuram(i):
    return int(subprocess.run(f"nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i {i}",
                 shell=True,capture_output=True,text=True).stdout.strip() or 99999)

log=open(f"{PBE}/run_all_pbe.log","a")
log.write(f"\n=== START @ {time.strftime('%H:%M:%S')} — {len(all_ins)} genuine-PBE calcs pending (priority: Y-decomp) ===\n"); log.flush()
print(f"{len(all_ins)} genuine-PBE calcs pending (priority: Y-decomp first)",flush=True)
running={}; done=0; total=len(all_ins)

def try_launch(gpu):
    if gpu in running or not all_ins: return
    if gpuram(gpu) > 2000: return   # GPU busy (another calc / convergence still running)
    inp=all_ins.pop(0); d=os.path.dirname(inp); tag=os.path.basename(inp)[:-3]
    for tmp in glob.glob(f"{d}/tmp_{tag}*"): subprocess.run(f"rm -rf {tmp}",shell=True)
    cmd=f"cd {d} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp {tag}.in > {tag}.out 2>&1"
    running[gpu]=(tag,subprocess.Popen(cmd,shell=True,env=env),time.time(),inp)
    msg=f"[{time.strftime('%H:%M:%S')}] GPU{gpu} START {tag.replace('_','/')}"; print(msg,flush=True); log.write(msg+"\n"); log.flush()

while running or all_ins:
    for gpu in (0,1): try_launch(gpu)
    time.sleep(15)
    for gpu in list(running):
        tag,p,t0,inp=running[gpu]
        if p.poll() is not None:
            dt=time.time()-t0; out=inp[:-3]+".out"
            txt=open(out,errors="ignore").read()
            ok="convergence has been achieved" in txt; pbe="Exchange-correlation= PBE" in txt
            done+=1
            msg=f"[{time.strftime('%H:%M:%S')}] GPU{gpu} DONE {tag.replace('_','/')} rc={p.returncode} conv={ok} pbe={pbe} ({dt:.0f}s) [{done}/{total}]"
            print(msg,flush=True); log.write(msg+"\n"); log.flush()
            for tmp in glob.glob(f"{os.path.dirname(inp)}/tmp_{tag}*"): subprocess.run(f"rm -rf {tmp}",shell=True)
            del running[gpu]
print(f"=== ALL {done} genuine-PBE SCFs COMPLETE ===",flush=True)
log.write(f"=== ALL {done} COMPLETE @ {time.strftime('%H:%M:%S')} ===\n"); log.flush()
