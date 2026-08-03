#!/usr/bin/env bash
# Inject HSE converged charge-density (+paw) into the complete PBE save of the SAME
# geometry, then pp.x plot_num=17 (PAW-reconstructed valence) and plot_num=0 (smooth
# valence). Sanity first: plot_num=17 on the PBE save must reproduce rho17.cube.
set -u
SP=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/hse800_spotcheck
YD=/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep
MPI=/workspace/openmpi; NVHPC=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11; CUDA=13.0
export PATH=$MPI/bin:$PATH
export LD_LIBRARY_PATH=$MPI/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/$CUDA/lib64:$NVHPC/cuda/$CUDA/lib64
QE=/workspace/qe-7.5
cd "$SP"
pp () {  # $1=plot_num $2=outfile $3=prefix $4=outdir
  cat > pp_tmp.in <<EOF
&INPUTPP
  prefix = '$3'
  outdir = '$4'
  plot_num = $1
  filplot = 'fld'
/
&PLOT
  nfile = 1
  filepp(1) = 'fld'
  weight(1) = 1.0
  iflag = 3
  output_format = 6
  fileout = '$2'
/
EOF
  CUDA_VISIBLE_DEVICES=0 $QE/bin/pp.x -inp pp_tmp.in > pp_$2.log 2>&1
  echo "  plot_num=$1 -> $2  exit=$?  $(ls -la $2 2>/dev/null | awk '{print $5"B"}')"
}
# --- template save (complete PBE, same geometry) ---
rm -rf tmpl_loaded tmpl_empty; mkdir -p tmpl_loaded tmpl_empty
cp -r "$YD/tmp_loaded/Y_deep_loaded.save" tmpl_loaded/
cp -r "$YD/tmp_empty/Y_deep_empty.save" tmpl_empty/ 2>/dev/null || echo "  (no PBE empty save; will handle CT separately)"
echo "=== template saves ready ==="; ls tmpl_loaded/ tmpl_empty/ 2>/dev/null
# --- SANITY: plot_num=17 on PBE density reproduces rho17.cube? ---
echo "=== SANITY: pp.x plot_num=17 on PBE save ==="
pp 17 sanity_rho17.cube Y_deep_loaded ./tmpl_loaded
# --- INJECT HSE density (+paw) into loaded template ---
echo "=== injecting HSE charge-density + paw into loaded template ==="
cp captured/hse800_vc.save/charge-density.dat tmpl_loaded/Y_deep_loaded.save/charge-density.dat
cp captured/hse800_vc.save/paw.txt          tmpl_loaded/Y_deep_loaded.save/paw.txt
echo "=== pp.x plot_num=17 (HSE) ==="
pp 17 hse_rho17.cube Y_deep_loaded ./tmpl_loaded
echo "=== pp.x plot_num=0 (HSE smooth valence, wfc1-independent) ==="
pp 0  hse_val.cube   Y_deep_loaded ./tmpl_loaded
echo "=== grids ==="
for f in sanity_rho17.cube hse_rho17.cube hse_val.cube; do [ -f "$f" ] && { printf "%-22s " "$f"; sed -n '4,6p' "$f" | awk '{printf "%s ",$1}'; echo; }; done
