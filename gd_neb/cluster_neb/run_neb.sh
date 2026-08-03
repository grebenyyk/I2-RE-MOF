#!/bin/bash
cd /workspace/dimitry/calc/qe/i2_mof/gd_neb/cluster_neb
echo "NEB START $(date)"
/workspace/miniconda3/bin/python3 run_neb_parallel.py 2>&1 | tee neb_run.log
echo "NEB END $(date)"
