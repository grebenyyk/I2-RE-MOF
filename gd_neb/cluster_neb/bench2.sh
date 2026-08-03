#!/bin/bash
source /workspace/miniconda3/etc/profile.d/conda.sh
conda activate nwchem
export NWCHEM_BASIS_LIBRARY=$CONDA_PREFIX/share/nwchem/libraries/
export NWCHEM_NWPW_LIBRARY=$CONDA_PREFIX/share/nwchem/libraryps/
export OMP_NUM_THREADS=1 OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
cd /workspace/dimitry/calc/qe/i2_mof/gd_neb/cluster_neb
T0=$(date +%s)
mpirun --allow-run-as-root --oversubscribe -np 16 nwchem sp_mol_pbe_svp.nw > sp_mol_pbe_svp.out 2>&1
T1=$(date +%s)
echo "WALL=$((T1-T0))s exit=$?"
