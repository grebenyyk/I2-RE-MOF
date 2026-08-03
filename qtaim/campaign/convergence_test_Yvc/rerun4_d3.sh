export PATH=/workspace/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/13.0/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2
QE=/workspace/qe-7.5/bin/pw.x
MISS="hmol_Yst_S1 hse_lmol_Yact_S2 hse_lmol_Yact_S3 hse_lpoly_Yst_S6"
i=0
for t in $MISS; do
  gpu=$([ $((i%2)) -eq 0 ] && echo 0 || echo 1); i=$((i+1))
  rm -rf tmp_d3only_$t
  CUDA_VISIBLE_DEVICES=$gpu mpirun --allow-run-as-root -np 1 "$QE" -nk 1 -nd 1 -nb 1 -nt 1 -inp d3only_$t.in > d3only_$t.out 2>&1 &
done
wait
echo "=== D3 results ==="
for t in $MISS; do
  d3=$(grep -E "DFT-D3 Dispersion +=" d3only_$t.out 2>/dev/null | tail -1 | awk '{print $NF}')
  echo "$t: D3=${d3:-MISSING}"
done
echo "RERUN4 DONE"
