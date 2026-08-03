#!/usr/bin/env python3
"""Wait for the 2 convergence calcs to finish (freeing GPUs), then relaunch
the full 69-calc genuine-PBE batch driver. Chained so GPUs are never contended."""
import subprocess,time,re,os
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/pbe_conv"
PBE="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results"
def done(tag):
    f=f"{D}/{tag}.out"
    if not os.path.exists(f): return False
    return "convergence has been achieved" in open(f,errors="ignore").read()
# wait for both convergence calcs
while not (done("pbe_k663") and done("pbe_poly_ec60")):
    time.sleep(30)
print(f"[{time.strftime('%H:%M:%S')}] convergence 9/9 done; waiting 15s for GPU release",flush=True)
time.sleep(15)
# relaunch batch driver
log=open(f"{PBE}/run_all_pbe.log","a")
log.write(f"\n=== RESTART after convergence @ {time.strftime('%H:%M:%S')} ===\n"); log.close()
subprocess.Popen(f"cd {PBE} && nohup python3 run_all_pbe.py >> run_all_pbe.log 2>&1 &",shell=True)
print(f"[{time.strftime('%H:%M:%S')}] batch driver restarted",flush=True)
