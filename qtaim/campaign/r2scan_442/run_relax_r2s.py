#!/usr/bin/env python3
"""Resumable driver: r2SCAN-NC I-only relax of Yact + Yst (2 systems, 2 GPUs).
- Idempotent: state from .out (converged marker) + <outdir>/<prefix>.save (resume point).
- Resume across node death: relaunch with restart_mode='restart' -> QE resumes relax
  from last ionic step (save dir persists on shared /workspace FS).
- GPU launch (CUDA_VISIBLE_DEVICES, NP=1, -nk1 -nd1 -nb1 -nt1), memory-guarded.
- Divergence guard: SCF acc >1 Ry @ iter 50 -> kill, flag DIVERGENT.
- Floor-acceptance: damp may not declare convergence on soft I modes; accept at
  force<0.04 + last-3 |dE|<0.06 kJ/mol (search-quality; final SP gives the number).
Usage: python3 run_relax_r2s.py <ecutwfc> <ecutrho>"""
import subprocess, time, os, re, fcntl, sys, glob
E  = sys.argv[1] if len(sys.argv) > 1 else "70"
RHO = sys.argv[2] if len(sys.argv) > 2 else str(int(E)*8)
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5-libxc-gpu"
env=os.environ.copy()
env["PATH"]=f"{MPI}/bin:{NVHPC}/compilers/bin:{NVHPC}/cuda/bin:{env.get('PATH','')}"
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/cuda/lib64:{env.get('LD_LIBRARY_PATH','')}"
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/r2scan_442"
TAGS=[f"rel_Yact_r2s_e{E}", f"rel_Yst_r2s_e{E}"]
MAXRETRY=6
LOG=open(f"{T}/run_relax.log","a")
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); LOG.write(line+"\n"); LOG.flush()
lf=open(f"{T}/relax.lock","w")
try: fcntl.flock(lf,fcntl.LOCK_EX|fcntl.LOCK_NB)
except (BlockingIOError,OSError): log("another driver instance holds the lock -- exiting"); sys.exit(0)

def gpuram(i): return int(subprocess.run(f"nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i {i}",shell=True,capture_output=True,text=True).stdout.strip() or 99999)
def outf(tag): return f"{T}/{tag}.out"
def inpf(tag): return f"{T}/{tag}.in"
def readf(p):
    try: return open(p,errors="ignore").read()
    except: return ""
def converged(tag):
    s=readf(outf(tag))
    if "bfgs converged" in s or "Damped Dynamics: convergence achieved" in s: return True
    # floor-acceptance (mirror search driver): >=15 ionic steps, last-3 |dE|<0.06 kJ, force<0.04
    if s.count("ATOMIC_POSITIONS")>=15:
        es=re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)",s)
        if len(es)>=4:
            dE=[abs(float(es[i])-float(es[i+1]))*13.6057 for i in range(len(es)-3,len(es)-1)]
            mf=re.findall(r"Total force\s*=\s*([0-9.]+)",s)
            if all(d<0.06 for d in dE) and mf and float(mf[-1])<0.04: return True
    return False
def failed(tag): return os.path.exists(f"{T}/{tag}.FAILED")
def retries(tag):
    f=f"{T}/{tag}.retries"; return int(open(f).read().strip()) if os.path.exists(f) else 0
def bump_retries(tag):
    f=f"{T}/{tag}.retries"; n=retries(tag)+1; open(f,"w").write(str(n)); return n
def save_exists(tag):
    base=open(inpf(tag)).read(); m=re.search(r"prefix\s*=\s*'([^']*)'",base); m2=re.search(r"outdir\s*=\s*'([^']*)'",base)
    return m and m2 and os.path.exists(f"{m2.group(1)}/{m.group(1)}.save")
def build_run_in(tag):
    base=open(inpf(tag)).read()
    # resume if a save exists (interrupted calc), else from_scratch. Keep outdir/prefix STABLE.
    if save_exists(tag) and "restart_mode" not in base:
        base=base.replace("&CONTROL","&CONTROL\n  restart_mode = 'restart'",1)
    open(f"{T}/{tag}_run.in","w").write(base); return f"{T}/{tag}_run.in"
def launch(gpu,tag):
    runf=build_run_in(tag)
    # do NOT rm the outdir -- the <prefix>.save in it is the resume point
    cmd=f"cd {T} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp {tag}_run.in > {tag}.out 2>&1"
    return subprocess.Popen(cmd,shell=True,env=env,start_new_session=True)
def scf_state(tag):
    raw=readf(outf(tag)); accs=re.findall(r"estimated scf accuracy\s*<\s*([0-9.]+)",raw); iters=re.findall(r"iteration #\s*(\d+)",raw)
    return (float(accs[-1]) if accs else None, int(iters[-1]) if iters else None, "ATOMIC_POSITIONS" in raw)

log(f"== driver start (ecutwfc={E} ecutrho={RHO}): {sum(converged(t) for t in TAGS)}/{len(TAGS)} converged ==")
running={}
while True:
    done=sum(converged(t) for t in TAGS); fail=sum(failed(t) for t in TAGS)
    if done+fail>=len(TAGS): log(f"== ALL RESOLVED: {done} converged, {fail} failed =="); break
    for g in (0,1):
        if g in running or gpuram(g)>3000: continue
        cand=[t for t in TAGS if not converged(t) and not failed(t) and t not in {running[gg]["tag"] for gg in running}]
        if not cand: continue
        tag=cand[0]; p=launch(g,tag); running[g]={"tag":tag,"p":p,"t0":time.time()}
        log(f"GPU{g} START {tag} {'resume' if save_exists(tag) else 'seed'} [{done}/{len(TAGS)}]")
    time.sleep(30)
    for g in list(running):
        e=running[g]; tag=e["tag"]; p=e["p"]
        if p.poll() is not None:
            if converged(tag):
                m=re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)",readf(outf(tag)))
                log(f"GPU{g} DONE {tag} E={m[-1] if m else '?'} ({time.time()-e['t0']:.0f}s)")
            else:
                acc,_,_=scf_state(tag); raw=readf(outf(tag))
                scf_failed=("convergence NOT been achieved" in raw) or ("Maximum CPU time" in raw)
                if scf_failed:
                    n=bump_retries(tag)
                    if n>=MAXRETRY: open(f"{T}/{tag}.FAILED","w").close(); log(f"GPU{g} {tag} SCF failed {n}x -> FAILED")
                    else: log(f"GPU{g} {tag} SCF failed -> retry {n}/{MAXRETRY} (will resume)")
                else: log(f"GPU{g} {tag} ended not-converged (nstep?) -> will resume")
            del running[g]
        else:
            if converged(tag):
                try: os.killpg(os.getpgid(p.pid),9)
                except: pass
                m=re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)",readf(outf(tag)))
                mf=re.findall(r"Total force\s*=\s*([0-9.]+)",readf(outf(tag)))
                log(f"GPU{g} DONE {tag} FLOOR-ACCEPT F={mf[-1][:6] if mf else '?'} E={m[-1] if m else '?'}"); del running[g]
            else:
                acc,it,_=scf_state(tag)
                if acc is not None and it is not None and acc>1.0 and it>50:
                    try: os.killpg(os.getpgid(p.pid),9)
                    except: pass
                    open(f"{T}/{tag}.FAILED","w").close(); log(f"GPU{g} {tag} KILL divergent >1Ry@it{it} -> FAILED"); del running[g]
log("driver exit")
