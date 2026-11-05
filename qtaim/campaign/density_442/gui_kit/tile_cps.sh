#!/usr/bin/env bash
# Locally tile a kit cube + regenerate its CP set to match (chains across cells).
# usage: tile_cps.sh <st|vc|act|actmol>   (needs python3 + numpy)
# st/vc/act -> [ -2,+2 ]x[ -1,+1 ]x[ -1,+1 ]; actmol -> [ -1,+1 ]^3
# NOTE: after this tiling the +a1 image of I3 (st plane) = atom 155.
set -euo pipefail
S=${1:?usage: tile_cps.sh <st|vc|act|actmol>}
case $S in
  st)     SRC=val_st.cube;     S0="0,1,-1,-2";;
  vc)     SRC=val_vc.cube;     S0="0,1,-1,-2";;
  act)    SRC=val_act.cube;    S0="0,1,-1,-2";;
  actmol) SRC=val_actmol.cube; S0="0,-1";;
esac
if [ "$S" = "actmol" ]; then
  python3 tile_cube.py "$SRC" t1.cube --axis 0 --neg
  python3 tile_cube.py t1.cube t2.cube --axis 1 --neg
  python3 tile_cube.py t2.cube val_${S}_big.cube --axis 2 --neg && rm t1.cube t2.cube
else
  python3 tile_cube.py "$SRC" t1.cube 2 --axis 0
  python3 tile_cube.py t1.cube t2.cube --axis 0 --neg
  python3 tile_cube.py t2.cube t1.cube --axis 1 --neg
  python3 tile_cube.py t1.cube val_${S}_big.cube --axis 2 --neg && rm t1.cube t2.cube
fi
python3 mwfn_cps_442.py rho17_${S}.cube CPs_${S}_big.txt \
  bcp_raw_${S}.out refine_${S}.out --shift0 $S0 --shift1 0,-1 --shift2 0,-1 \
  2>/dev/null || echo "(rho17/bcp_raw/refine not shipped in the kit - reuse unit CPs.txt, \
dots outside cell 0 will be missing: cp CPs_${S}.txt CPs_${S}_big.txt)"
echo "wrote val_${S}_big.cube (atoms 1-80 unchanged; +a1 image block starts at 81)"
