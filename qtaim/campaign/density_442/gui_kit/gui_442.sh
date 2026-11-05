#!/usr/bin/env bash
# Launch Multiwfn in GUI mode (interactive graph windows) on THIS node via X11
# forwarding: from your desktop connect with  ssh -X  (or MobaXterm/VcXsrv),
# then:  ./gui_442.sh val_st.cube CPs_st.txt
# Each plane map pops a native window; close it to continue in the menu.
set -uo pipefail
BIN=/workspace/dimitry/calc/qe/i2_mof/tools/multiwfn/Multiwfn_2026.7.11_bin_Linux/Multiwfn
export MultiwfnPATH="$(dirname $BIN)" OMP_STACKSIZE=128m
CUBE=${1:?usage: gui_442.sh <cube> [CPs.txt]}; CPS=${2:-CPs.txt}
if [ -z "${DISPLAY:-}" ]; then echo "No DISPLAY - connect with ssh -X (or set DISPLAY=:0)"; exit 1; fi
rm -f CPs.txt settings.ini
cp "$CPS" CPs.txt
sed -e 's/^  isilent=.*/  isilent= 0/' -e 's/^  graphformat=.*/  graphformat= png/' \
  /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/bcp_plots/multiwfn_out/settings.ini > settings.ini
"$BIN"
