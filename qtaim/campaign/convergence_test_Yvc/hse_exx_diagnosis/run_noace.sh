#!/usr/bin/env bash
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/hse_exx_diagnosis
export PATH=/workspace/miniconda3/envs/qe/bin:$PATH
export LD_LIBRARY_PATH=/workspace/miniconda3/envs/qe/lib:${LD_LIBRARY_PATH:-}
rm -rf tmp
mpirun --allow-run-as-root -np 1 pw.x -inp gas_hse_noace.in > gas_hse_noace.out 2>&1
echo "EXIT=$?" >> gas_hse_noace.out
