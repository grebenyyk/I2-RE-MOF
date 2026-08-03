#!/bin/bash
set -u
D=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc
QE=/workspace/qe-7.5/bin/pw.x
export PATH=/workspace/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/workspace/openmpi/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/math_libs/13.0/lib64:/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/cuda/13.0/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2
echo $$ > "$D/_run/relax442.pid"
echo "=== relax442 start $(date) ===" >> "$D/_run/relax442_status.log"
# GPU0: molecular start
( cd "$D" && rm -rf tmp_relax442_mol_yact && CUDA_VISIBLE_DEVICES=0 mpirun --allow-run-as-root -np 1 "$QE" -nk 1 -nd 1 -nb 1 -nt 1 -inp relax442_mol_yact.in > relax442_mol_yact.out 2>&1 && touch _run/relax442_mol_yact.done || touch _run/relax442_mol_yact.fail ) &
# GPU1: polyiodide start
( cd "$D" && rm -rf tmp_relax442_poly_yact && CUDA_VISIBLE_DEVICES=1 mpirun --allow-run-as-root -np 1 "$QE" -nk 1 -nd 1 -nb 1 -nt 1 -inp relax442_poly_yact.in > relax442_poly_yact.out 2>&1 && touch _run/relax442_poly_yact.done || touch _run/relax442_poly_yact.fail ) &
wait
echo "=== relax442 finished $(date) ===" >> "$D/_run/relax442_status.log"
