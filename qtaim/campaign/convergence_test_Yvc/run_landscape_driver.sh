#!/bin/bash
# Landscape pipeline driver. Called by cron. Handles:
# 1. GPU1 batch launch (when mmol_Yst finishes)
# 2. Plateau detection (mark energy-converged relaxes as done)
# 3. HSE build+launch (when all relaxed)
# 4. Analysis (when all HSE done)
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc
QE=/workspace/qe-7.5/bin/pw.x
export PATH=/workspace/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/13.0/lib64:${LD_LIBRARY_PATH:-}
cd "$D"
ALL="mmol_Yst lmol_Yact_S2 lmol_Yact_S3 lmol_Yvc_S2 lmol_Yvc_S3 lmol_Yst_S2 lmol_Yst_S3 lpoly_Yact_S1 lpoly_Yact_S2 lpoly_Yact_S3 lpoly_Yvc_S2 lpoly_Yvc_S3 lpoly_Yst_S2 lpoly_Yst_S3"
# 1. GPU1 batch launch (when mmol_Yst done)
if [ ! -f _run/landscape_gpu1_launched.marker ] && [ -f _run/mmol_Yst.done ]; then
  touch _run/landscape_gpu1_launched.marker
  nohup bash run_mseeds.sh 1 lmol_Yvc_S3 lmol_Yst_S2 lpoly_Yvc_S3 lmol_Yst_S3 lpoly_Yst_S2 lpoly_Yst_S3 >run_landscape_gpu1.log 2>&1 &
  echo "[$(date +%H:%M)] GPU1 batch launched (6 relaxes)"
fi
# 2. plateau detection: mark energy-converged relaxes as done
for tag in $ALL; do
  [ -f "_run/${tag}.done" ] && continue
  [ ! -f "${tag}.out" ] && continue
  steps=$(grep -c 'Forces acting on atoms' "${tag}.out" 2>/dev/null)
  force=$(grep 'Total force' "${tag}.out" 2>/dev/null | tail -1 | awk '{print $4}')
  # check: relax process not running (pw.x not on this tag) + force<0.05 + steps>=8
  if [ -n "$force" ] && python3 -c "exit(0 if (float('$force')<0.05 and $steps>=8) or (float('$force')<0.2 and $steps>=10) else 1)" 2>/dev/null; then
    # plateau: kill pw.x to stop wasting GPU time, then mark done
    for p in $(pgrep -x pw.x 2>/dev/null); do
      cl=$(tr '\0' ' ' </proc/$p/cmdline 2>/dev/null)
      case "$cl" in *${tag}*) kill "$p" 2>/dev/null;; esac
    done
    sleep 2
    touch "_run/${tag}.done"
    echo "[$(date +%H:%M)] $tag plateau-killed+done (force=$force, $steps steps)"
  fi
done
# 3. count relaxed + report
nrelaxed=0; ntotal=14
for tag in $ALL; do [ -f "_run/${tag}.done" ] && nrelaxed=$((nrelaxed+1)); done
echo "relaxed: $nrelaxed/$ntotal"
for tag in $ALL; do
  d=$([ -f "_run/${tag}.done" ] && echo "✓" || echo ".")
  h=$([ -f "_run/hse_${tag}.done" ] && echo "✓" || echo ".")
  printf "  %-16s relax=%s hse=%s\n" "$tag" "$d" "$h"
done
# 4. when all relaxed, build+launch HSE (once)
if [ "$nrelaxed" = "$ntotal" ] && [ ! -f _run/hse_all_launched.marker ]; then
  touch _run/hse_all_launched.marker
  echo "[$(date +%H:%M)] ALL RELAXED. Building+launching HSE..."
  hse_list=""
  for tag in $ALL; do
    hse_tag="hse_${tag}"
    python3 build_mol_hse_single.py "$tag" "$hse_tag" 2>/dev/null
    hse_list="$hse_list $hse_tag"
  done
  # split HSE across 2 GPUs
  i=0
  for h in $hse_list; do
    gpu=$([ $((i%2)) -eq 0 ] && echo 0 || echo 1); i=$((i+1))
    nohup bash run_hse_conv.sh "$gpu" "$h" >"run_${h}.log" 2>&1 &
    sleep 2
  done
  echo "[$(date +%H:%M)] HSE launched ($ntotal calcs, 2 GPUs)"
fi
# 4b. re-launch any HSE calcs that failed/OOM'd (when GPU memory frees up)
for tag in $ALL; do
  hse_tag="hse_${tag}"
  if [ ! -f "_run/${hse_tag}.done" ] && [ ! -f "_run/${hse_tag}_relaunched.marker" ]; then
    running=$(for p in $(pgrep -x pw.x 2>/dev/null); do tr '\0' ' ' </proc/$p/cmdline 2>/dev/null; done | grep -c "$hse_tag")
    if [ "$running" = "0" ]; then
      g0=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | head -1)
      g1=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | tail -1)
      if [ -n "$g0" ] && [ "$g0" -lt 25000 ]; then gpu=0
      elif [ -n "$g1" ] && [ "$g1" -lt 25000 ]; then gpu=1
      else break; fi
      touch "_run/${hse_tag}_relaunched.marker"
      nohup bash run_hse_conv.sh "$gpu" "$hse_tag" >"run_${hse_tag}_rerun.log" 2>&1 &
      sleep 5
      echo "[$(date +%H:%M)] re-launched $hse_tag on GPU$gpu (mem g0=$g0 g1=$g1)"
    fi
  fi
done
# 5. when all HSE done, analyze
nhse=0
for tag in $ALL; do [ -f "_run/hse_${tag}.done" ] && nhse=$((nhse+1)); done
if [ "$nhse" = "$ntotal" ] && [ ! -f _run/landscape_analyzed.marker ]; then
  touch _run/landscape_analyzed.marker
  echo "[$(date +%H:%M)] ALL HSE DONE. Running analysis..."
  python3 analyze_landscape_full.py
fi
