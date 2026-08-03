#!/usr/bin/env python3
"""Run all campaign SCFs (loaded+empty+i2only per system) on 2 GPUs. Restartable (.done markers)."""
import subprocess, os, time
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
QE="/workspace/qe-7.5"; MPI="/workspace/openmpi"
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; CUDA="13.0"
os.environ["PATH"]=MPI+"/bin:"+os.environ["PATH"]
os.environ["LD_LIBRARY_PATH"]=MPI+"/lib:"+NVHPC+"/compilers/lib:"+NVHPC+f"/math_libs/{CUDA}/lib64:"+NVHPC+f"/cuda/{CUDA}/lib64"
SYSTEMS=["Gd_deep","Tb_deep","Y_deep"]
CALCS=["loaded","empty","i2only"]
def marker(s,c): return f"{BASE}/{s}/{c}.done"
def conv(s,c):
    p=f"{BASE}/{s}/scf_{c}.out"
    return "convergence has been achieved" in open(p,errors="ignore").read() if os.path.exists(p) else False
jobs=[(s,c) for s in SYSTEMS for c in CALCS]
jobs=[j for j in jobs if not (os.path.exists(marker(*j)) or conv(*j))]
print(f"{len(jobs)} SCFs to run (of {len(SYSTEMS)*len(CALCS)})",flush=True)
running={}  # gpu -> (sys,calc,popen,cmdstr)
def launch(gpu):
    s,c=jobs.pop(0); d=f"{BASE}/{s}"
    cmd=f"cd {d} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp scf_{c}.in > scf_{c}.out 2>&1"
    running[gpu]=(s,c,subprocess.Popen(cmd,shell=True),time.time())
    print(f"[{time.strftime('%H:%M:%S')}] GPU{gpu} START {s}/{c}",flush=True)
log=open(f"{BASE}/scf_campaign.log","a")
while jobs or running:
    for gpu in (0,1):
        if gpu not in running and jobs:
            launch(gpu)
    time.sleep(30)
    for gpu in list(running):
        s,c,p,t0=running[gpu]
        if p.poll() is not None:
            ok=conv(s,c)
            if ok: open(marker(s,c),"w").write("ok")
            dt=time.time()-t0
            msg=f"[{time.strftime('%H:%M:%S')}] GPU{gpu} DONE {s}/{c} rc={p.returncode} conv={ok} ({dt:.0f}s)"
            print(msg,flush=True); log.write(msg+"\n"); log.flush()
            del running[gpu]
print("=== ALL SCFs COMPLETE ===",flush=True)
# summary
for s in SYSTEMS:
    for c in CALCS:
        print(f"  {s}/{c}: {'DONE' if os.path.exists(marker(s,c)) else ('CONV' if conv(s,c) else 'MISSING')}",flush=True)
