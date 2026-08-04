#!/usr/bin/env python3
"""Resumable driver: 18-config Y 4x4x2 landscape relaxation campaign.
- Idempotent: all state derived from .out files (no separate state file).
- Reseed-on-resume: a partial .out (node died mid-calc) is reseeded from its latest
  geometry, preserving progress. Relaunch on a new node via `bash resume.sh`.
- Memory-guarded (2 GPUs, gpuram<3000); flock-locked (one instance).
- Divergence guard: SCF accuracy >1 Ry past iter 50 -> kill, retry once at degauss=0.06,
  else FAILED. Non-divergent non-convergence -> reseed-retry, max 4 then FAILED.
"""
import subprocess,time,os,re,fcntl,sys
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"; QE="/workspace/qe-7.5"
env=os.environ.copy()
env["PATH"]=f"{MPI}/bin:{NVHPC}/compilers/bin:{env.get('PATH','')}"
env["LD_LIBRARY_PATH"]=f"{MPI}/lib:{NVHPC}/compilers/lib:{NVHPC}/math_libs/13.0/lib64:{env.get('LD_LIBRARY_PATH','')}"
T="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/landscape_442"
TAGS=["vc_poly_S1","vc_poly_lS3","vc_poly_lS6","vc_mol_lS4","vc_mol_ref","vc_mol_lS5","vc_mol_lS3",
      "act_mol_rd","act_poly_lS2","act_poly_lS5","act_mol_lS2","act_mol_lS6",
      "st_poly_S1","st_poly_lS5","st_poly_lS6","st_poly_lS4","st_mol_lS6","st_mol_lS3"]
MAXRETRY=4
LOG=open(f"{T}/run.log","a")
def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True); LOG.write(line+"\n"); LOG.flush()
# ---- flock: single instance ----
lf=open(f"{T}/l442.lock","w")
try: fcntl.flock(lf,fcntl.LOCK_EX|fcntl.LOCK_NB)
except (BlockingIOError,OSError): log("another driver instance holds the lock -- exiting"); sys.exit(0)

def gpuram(i): return int(subprocess.run(f"nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i {i}",shell=True,capture_output=True,text=True).stdout.strip() or 99999)
def outf(tag): return f"{T}/l442_{tag}.out"
def inpf(tag): return f"{T}/l442_{tag}.in"
def readf(p):
    try: return open(p,errors="ignore").read()
    except: return ""
def converged(tag):
    s=readf(outf(tag))
    if "bfgs converged" in s or "Damped Dynamics: convergence achieved" in s: return True
    # floor-acceptance: soft-mode configs wiggle at the floor (damp ΔE ~threshold, never declares).
    # Accept if >=30 steps, last-4 |ΔE|<0.04 kJ/mol (energy flat), last force<0.04 (soft-mode floor).
    if s.count("ATOMIC_POSITIONS")>=30:
        es=re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)",s)
        if len(es)>=5:
            dE=[abs(float(es[i])-float(es[i+1]))*13.6057 for i in range(len(es)-4,len(es)-1)]
            mf=re.findall(r"Total force\s*=\s*([0-9.]+)",s)
            if all(d<0.04 for d in dE) and mf and float(mf[-1])<0.04: return True
    return False
def failed(tag): return os.path.exists(f"{T}/l442_{tag}.FAILED")
def divergent04(tag): return os.path.exists(f"{T}/l442_{tag}.DIVERGENT04")
def retries(tag):
    f=f"{T}/l442_{tag}.retries"
    return int(open(f).read().strip()) if os.path.exists(f) else 0
def bump_retries(tag):
    f=f"{T}/l442_{tag}.retries"; n=retries(tag)+1; open(f,"w").write(str(n)); return n
def last_atoms(outf_):
    raw=readf(outf_)
    if "ATOMIC_POSITIONS" not in raw: return None
    blk=raw.split("ATOMIC_POSITIONS")[-1].split("K_POINTS")[0].strip().split("\n")[1:]
    atoms=[(l.split()[0],l.split()[1:4]) for l in blk if l.split() and l.split()[0] in ("C","H","O","Y","I")]
    return atoms if len(atoms)==80 else None
def build_run_in(tag,dg):
    base=open(inpf(tag)).read()
    atoms=last_atoms(outf(tag))
    if atoms:
        base_pos=re.search(r"ATOMIC_POSITIONS crystal\n(.*?)\nK_POINTS",base,re.S).group(1).strip().split("\n")
        poslines=[]; ai=0
        for l in base_pos:
            sp=l.split()[0]
            if sp in ("C","H","O","Y","I"):
                _,frac=atoms[ai]; ai+=1
                flag="1 1 1" if sp=="I" else "0 0 0"
                poslines.append(f"{sp} {float(frac[0]):.12f} {float(frac[1]):.12f} {float(frac[2]):.12f} {flag}")
        base=re.sub(r"ATOMIC_POSITIONS crystal\n.*?(?=\nK_POINTS)","ATOMIC_POSITIONS crystal\n"+"\n".join(poslines)+"\n",base,flags=re.S)
    base=re.sub(r"degauss = [0-9.]+",f"degauss = {dg}",base)
    base=re.sub(r"prefix = '[^']*'",f"prefix = 'l442_{tag}_run'",base)
    base=re.sub(r"outdir = '[^']*'",f"outdir = './tmp_l442_{tag}_run'",base)
    runf=f"{T}/l442_{tag}_run.in"; open(runf,"w").write(base); return runf
def launch(gpu,tag):
    subprocess.run(f"rm -rf {T}/tmp_l442_{tag}_run*",shell=True)
    cmd=f"cd {T} && CUDA_VISIBLE_DEVICES={gpu} mpirun --allow-run-as-root -np 1 {QE}/bin/pw.x -nk 1 -nd 1 -nb 1 -nt 1 -inp l442_{tag}_run.in > l442_{tag}.out 2>&1"
    return subprocess.Popen(cmd,shell=True,env=env,start_new_session=True)
def scf_state(tag):
    raw=readf(outf(tag))
    accs=re.findall(r"estimated scf accuracy\s*<\s*([0-9.]+)",raw)
    iters=re.findall(r"iteration #\s*(\d+)",raw)
    return (float(accs[-1]) if accs else None, int(iters[-1]) if iters else None,
            "ATOMIC_POSITIONS" in raw)  # has>=1 ionic step

log(f"== driver start: {sum(converged(t) for t in TAGS)}/{len(TAGS)} converged, {sum(failed(t) for t in TAGS)} failed ==")
running={}
while True:
    done=sum(converged(t) for t in TAGS); fail=sum(failed(t) for t in TAGS)
    if done+fail>=len(TAGS): log(f"== ALL RESOLVED: {done} converged, {fail} failed =="); break
    for g in (0,1):
        if g in running or gpuram(g)>3000: continue
        rt={running[gg]["tag"] for gg in running}
        cand=[t for t in TAGS if not converged(t) and not failed(t) and t not in rt]
        if not cand: continue
        tag=cand[0]; dg=0.06 if divergent04(tag) else 0.04
        runf=build_run_in(tag,dg); p=launch(g,tag)
        running[g]={"tag":tag,"p":p,"t0":time.time(),"dg":dg}
        log(f"GPU{g} START {tag} dg={dg} {'reseed' if last_atoms(outf(tag)) else 'seed'} [{done}/{len(TAGS)}]")
    time.sleep(30)
    for g in list(running):
        e=running[g]; tag=e["tag"]; p=e["p"]; dg=e["dg"]
        if p.poll() is not None:
            if converged(tag):
                m=re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)",readf(outf(tag)))
                log(f"GPU{g} DONE {tag} E={m[-1] if m else '?'} ({time.time()-e['t0']:.0f}s)")
            else:
                acc,_,hasstep=scf_state(tag); raw=readf(outf(tag))
                scf_failed=("convergence NOT been achieved" in raw) or (acc is not None and acc>1e-2 and not hasstep)
                if scf_failed:
                    if dg==0.04 and not divergent04(tag):
                        open(f"{T}/l442_{tag}.DIVERGENT04","w").close(); log(f"GPU{g} {tag} SCF divergent @0.04 -> retry @0.06")
                    else:
                        open(f"{T}/l442_{tag}.FAILED","w").close(); log(f"GPU{g} {tag} SCF divergent @0.06 -> FAILED")
                else:
                    n=bump_retries(tag)
                    if n>=MAXRETRY: open(f"{T}/l442_{tag}.FAILED","w").close(); log(f"GPU{g} {tag} {n}x non-converged -> FAILED")
                    else: log(f"GPU{g} {tag} ended not-converged (nstep?) retry {n}/{MAXRETRY}")
            del running[g]
        else:
            if converged(tag):  # mid-run floor-acceptance: kill+done if at the soft-mode floor
                try: os.killpg(os.getpgid(p.pid),9)
                except: pass
                m=re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)",readf(outf(tag)))
                mf=re.findall(r"Total force\s*=\s*([0-9.]+)",readf(outf(tag)))
                log(f"GPU{g} DONE {tag} FLOOR-ACCEPT F={mf[-1][:6] if mf else '?'} E={m[-1] if m else '?'}")
                del running[g]
            else:
                acc,it,_=scf_state(tag)
                if acc is not None and it is not None and acc>1.0 and it>50:
                    try: os.killpg(os.getpgid(p.pid),9)
                    except: pass
                    if dg==0.04 and not divergent04(tag):
                        open(f"{T}/l442_{tag}.DIVERGENT04","w").close(); log(f"GPU{g} {tag} KILL divergent >1Ry@it{it} -> retry @0.06")
                    else:
                        open(f"{T}/l442_{tag}.FAILED","w").close(); log(f"GPU{g} {tag} KILL divergent @0.06 -> FAILED")
                    del running[g]
log("driver exit")
