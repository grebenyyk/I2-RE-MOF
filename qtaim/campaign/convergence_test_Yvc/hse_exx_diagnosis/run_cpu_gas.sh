#!/usr/bin/env bash
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc/hse_exx_diagnosis
export PATH=/workspace/openmpi/bin:/workspace/miniconda3/envs/qe/bin:$PATH
export LD_LIBRARY_PATH=/workspace/miniconda3/envs/qe/lib:/workspace/openmpi/lib:${LD_LIBRARY_PATH:-}
rm -rf tmp
mpirun --allow-run-as-root -np 4 /workspace/miniconda3/envs/qe/bin/pw.x -inp gas_hse.in > gas_hse_cpu.out 2>&1
echo "EXIT=$?" >> gas_hse_cpu.out
