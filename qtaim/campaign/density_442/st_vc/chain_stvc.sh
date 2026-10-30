#!/bin/bash
# waits for SCF batch, then cubes + deltarho + critic2 for st/vc
T=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/st_vc
cd $T
while ! grep -q "rs_stvc batch END" run_rs_stvc.log 2>/dev/null; do sleep 120; done
echo "SCFs done, generating cubes $(date +%H:%M:%S)"
NVHPC="/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"; MPI="/workspace/openmpi"
export PATH="$MPI/bin:$NVHPC/compilers/bin:$PATH"
export LD_LIBRARY_PATH="$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/13.0/lib64:${LD_LIBRARY_PATH:-}"
for tag in st_poly_lS5 vc_poly_lS3; do
  for kind,pl,fn in dummy; do :; done 2>/dev/null
done
python3 <<'PYEOF'
import os
T='/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/st_vc/'
jobs=[]
for tag in ('st_poly_lS5','vc_poly_lS3'):
    pre={'loaded':'rs_loaded_','empty':'rs_empty_','i2only':'rs_i2only_'}
    for kind,plot,fn in [('loaded',17,f'rho17_{tag}.cube'),('loaded',0,f'val_{tag}.cube'),('loaded',8,f'elf_{tag}.cube'),
                          ('empty',0,f'empty_val_{tag}.cube'),('i2only',0,f'i2only_val_{tag}.cube')]:
        jobs.append((pre[kind]+tag,plot,fn))
for src,plot,fn in jobs:
    p=T+f"pp_{fn.replace('.cube','')}.in"
    open(p,'w').write(f"""&INPUTPP
  prefix = '{src}'
  outdir = './tmp_{src}'
  plot_num = {plot}
  filplot = 'f'
/
&PLOT
  nfile = 1
  filepp(1) = 'f'
  weight(1) = 1.0
  iflag = 3
  output_format = 6
  fileout = '{fn}'
/
""")
print('pp inputs:',len(jobs))
PYEOF
for f in pp_*.in; do c=$(grep fileout $f | awk '{print $3}' | tr -d "'"); [ -f "$c" ] && continue
  mpirun --allow-run-as-root -np 4 /workspace/qe-7.5/bin/pp.x -inp $f > ${f%.in}.out 2>&1 && echo "OK $c" || echo "FAIL $c"; done
python3 <<'PYEOF'
import numpy as np
T='/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/density_442/st_vc/'
def read_cube(p):
    L=open(p).readlines()
    nat=int(L[2].split()[0]); dims=[int(L[i].split()[0]) for i in (3,4,5)]
    data=np.loadtxt(p,skiprows=6+nat)
    return dims,data.reshape(dims)
def write_cube_like(src,data,out):
    L=open(src).readlines()[:6+int(open(src).readlines()[2].split()[0])]
    with open(out,'w') as g:
        g.writelines(L); np.savetxt(g,data.ravel(),fmt='%.5e')
for tag in ('st_poly_lS5','vc_poly_lS3'):
    d,vl=read_cube(T+f'val_{tag}.cube'); _,ve=read_cube(T+f'empty_val_{tag}.cube'); _,vi=read_cube(T+f'i2only_val_{tag}.cube')
    assert d[0]>=150
    write_cube_like(T+f'val_{tag}.cube',vl-ve-vi,T+f'deltarho_{tag}.cube')
    print(f'deltarho_{tag}: net={((vl-ve-vi).sum()):+.2e}')
    open(T+f'bcp_{tag}.cri','w').write(f"crystal rho17_{tag}.cube\nload rho17_{tag}.cube id rho\nauto\n")
    open(T+f'ct_{tag}.cri','w').write(f"crystal rho17_{tag}.cube\nload rho17_{tag}.cube id rho\nload deltarho_{tag}.cube id dr\nintegrable dr\nbader\n")
print('cri written')
PYEOF
export OMP_NUM_THREADS=1
C2=/workspace/dimitry/calc/qe/i2_mof/tools/critic2_run.sh
for j in bcp_st_poly_lS5 ct_st_poly_lS5 bcp_vc_poly_lS3 ct_vc_poly_lS3; do
  $C2 $j.cri > $j.out 2>&1 && echo "OK $j" || echo "FAIL $j"
done
echo "STVC CHAIN ALL DONE"
