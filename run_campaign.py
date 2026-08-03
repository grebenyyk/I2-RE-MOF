#!/usr/bin/env python3
"""
Run the I2@MOF campaign on 2 GPUs without oversubscription.

Strategy: 2 worker threads, each pinned to ONE GPU via CUDA_VISIBLE_DEVICES,
each running one pw.x at NP=1 (1 rank -> 1 GPU). So exactly 2 calcs run at a
time, one per GPU -- no VRAM oversubscription (safe even for the 206-atom Cd
cells). A shared queue keeps both GPUs fed: empties (fast SCF) first, then the
slow loaded relaxes.

Each job calls run.sh in the FOREGROUND (no --detach) so the worker blocks
until pw.x exits and knows the return code.

Run detached:  nohup python run_campaign.py > campaign.log 2>&1 &
Then poll:     python collect_energies.py ; tail campaign.log
"""
import os, subprocess, threading, time, queue

BASE = "/workspace/dimitry/calc/qe/i2_mof"
# MOFS and GPU ids overridable via env, e.g.:
#   MOFS="10_Cd 12_Tb 17_Cd" GPU_IDS="0,1" python run_campaign.py
_default_mofs = ["05_Cd", "06_Tb", "07_Cd", "10_Cd", "12_Tb", "17_Cd"]
MOFS = os.environ.get("MOFS", " ".join(_default_mofs)).split()
GPU_IDS = [int(x) for x in os.environ.get("GPU_IDS", "0,1").split(",") if x.strip()]
NGPU = len(GPU_IDS)

# If JOBS env is set (space-separated input paths, output = inp with .in->.out), use those
# directly instead of the MOFS empty+loaded expansion. E.g. for variant loaded-only runs.
if os.environ.get("JOBS"):
    _job_inps = os.environ["JOBS"].split()
    JOBS = [(os.path.dirname(inp) or ".", os.path.basename(inp).replace(".in", ""),
             inp, inp.replace(".in", ".out")) for inp in _job_inps]
else:
    JOBS = ([(m, "empty_mof", f"{m}/empty_mof.in", f"{m}/empty_mof.out") for m in MOFS] +
            [(m, "loaded",    f"{m}/loaded.in",    f"{m}/loaded.out")    for m in MOFS])

q = queue.Queue()
for j in JOBS:
    q.put(j)
lock = threading.Lock()
results = []


def log(msg):
    with lock:
        print(msg, flush=True)


def worker(gpu):
    env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(gpu), NP="1", OMP_NUM_THREADS="1")
    while True:
        try:
            m, kind, inp, out = q.get_nowait()
        except queue.Empty:
            return
        outp = os.path.join(BASE, out)
        if os.path.isfile(outp) and "JOB DONE" in open(outp, errors="ignore").read():
            log(f"[GPU{gpu}] SKIP  {m}/{kind} (already JOB DONE)")
            with lock:
                results.append((m, kind, 0, 0.0))
            continue
        t0 = time.time()
        log(f"[GPU{gpu}] START {m}/{kind}")
        rc = subprocess.run(["bash", "run.sh", inp, out], cwd=BASE, env=env).returncode
        dt = time.time() - t0
        with lock:
            results.append((m, kind, rc, dt))
        log(f"[GPU{gpu}] DONE  {m}/{kind}  rc={rc}  ({dt/60:.1f} min)")


def main():
    log(f"campaign: {len(JOBS)} jobs, {NGPU} GPUs, 1 calc/GPU (NP=1)")
    log("queue: " + ", ".join(f"{m}/{k}" for m, k, _, _ in JOBS))
    t0 = time.time()
    threads = [threading.Thread(target=worker, args=(g,)) for g in GPU_IDS]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    log(f"=== ALL DONE in {(time.time()-t0)/60:.1f} min ===")
    ok = sum(1 for _, _, rc, _ in results if rc == 0)
    log(f"summary: {ok}/{len(results)} jobs rc=0")
    for m, kind, rc, dt in results:
        log(f"  {m:6s}/{kind:9s} rc={rc} {dt/60:6.1f}min")


if __name__ == "__main__":
    main()
