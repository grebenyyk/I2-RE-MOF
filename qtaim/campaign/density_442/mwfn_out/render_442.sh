#!/usr/bin/env bash
# Render one 442 QTAIM plane map
# usage: render_442.sh <cube> <atoms> <CPs> <out.png|pdf> [extension_bohr]
# graphformat taken from the output extension. Post-menu kept MINIMAL:
# 8+3 bonds blue, -8 Angstrom, 0 save. Do NOT press 17/18 (label
# style/threshold): this build toggles labels off unpredictably; settings.ini
# defaults give element-symbol labels for atoms/CPs within 0.5 bohr of the
# plane. Harmless stdin forrtl note at "q".
set -uo pipefail
BIN=/workspace/dimitry/calc/qe/i2_mof/tools/multiwfn/Multiwfn_2026.7.11_bin_Linux/Multiwfn
export MultiwfnPATH="$(dirname $BIN)" OMP_STACKSIZE=128m
CUBE=$1; ATOMS=$2; CPS=$3; OUT=$4; EXT=${5:-10}
FMT=${OUT##*.}
rm -f CPs.txt dislin.png dislin.pdf settings.ini
sed "s/^  graphformat=.*/  graphformat= $FMT/" \
  /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/bcp_plots/multiwfn_out/settings.ini > settings.ini
cp "$CPS" CPs.txt
printf '%s\n\n2\n-4\n5\n\n0\n-10\n4\n100\n6\n300,300\n0\n%s\n4\n%s\n8\n3\n-8\n0\nq\n' \
  "$CUBE" "$EXT" "$ATOMS" | timeout 900 xvfb-run -a "$BIN" > "log_$(basename $OUT .${FMT}).out" 2>&1 || true
test -f "dislin.$FMT" || { echo "RENDER FAILED - see log_$(basename $OUT .${FMT}).out"; exit 1; }
mv "dislin.$FMT" "$OUT"
echo "wrote $OUT"
