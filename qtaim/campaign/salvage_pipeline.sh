#!/bin/bash
# Full HSE-PDOS salvage pipeline for one variant:
#   HSE maxstep=30 (get .wfc1 + schema) -> PBE scf (template wfc1.dat) -> decode -> projwfc
# usage: salvage_pipeline.sh <gpu> <variant_dir>
set -e
export MPI=/workspace/openmpi
export NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11
export CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
export OMP_NUM_THREADS=4
GPU=$1; DIR=$2; name=$(basename $DIR); prefix=${name}_hse
QE=/workspace/qe-7.5/bin/pw.x
CAMP=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign
cd "$DIR"
echo "[$name] step 1: HSE maxstep=30 -> .wfc1 + schema"
rm -rf tmp_hse
CUDA_VISIBLE_DEVICES=$GPU mpirun --allow-run-as-root -np 1 $QE -nk 1 -nd 1 -nb 1 -nt 1 -inp hse/scf.in > hse/scf.out 2>&1
mkdir -p hse/hse_keep; rm -rf hse/hse_keep/save
cp -r tmp_hse/${prefix}.save hse/hse_keep/save
cp tmp_hse/${prefix}.wfc1 hse/hse_keep/hse_local.wfc1
echo "[$name] step 2: PBE scf -> template"
rm -rf tmp_hse
CUDA_VISIBLE_DEVICES=$GPU mpirun --allow-run-as-root -np 1 $QE -nk 1 -nd 1 -nb 1 -nt 1 -inp hse/scf_pbe.in > hse/scf_pbe.out 2>&1
rm -rf hse/pbe_template_save; cp -r tmp_hse/${prefix}.save hse/pbe_template_save
echo "[$name] step 3: decode"
cd "$CAMP"; python3 decode_hse_pdos.py "$DIR"
echo "[$name] step 4: projwfc"
cd "$DIR"; $QE/../../bin/projwfc.x -inp hse/proj_decoded.in > hse/proj_decoded.out 2>&1 || /workspace/qe-7.5/bin/projwfc.x -inp hse/proj_decoded.in > hse/proj_decoded.out 2>&1
echo "[$name] DONE: $(ls hse_pdos* 2>/dev/null | wc -l) pdos files"
