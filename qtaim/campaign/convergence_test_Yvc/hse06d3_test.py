#!/usr/bin/env python3
"""Convert a d3only_*.in (PBE+D3) into an HSE06-D3 variant: input_dft='hse' +
exx_fraction=0 so the SCF is PBE-cost but the D3 functional is hse06 (QE maps
hse->hse06). D3 is geometry-only so ecut is irrelevant. Usage: conv.py <in> <out>"""
import sys,re
src=open(sys.argv[1]).read()
# inject hse+exx0 cards into &SYSTEM (after dftd3_version line)
add="""  input_dft='hse'
  exx_fraction=0.0d0
  screening_parameter=0.106d0
  exxdiv_treatment='gygi-baldereschi'
  x_gamma_extrapolation=.true.
  ecutfock=100
"""
assert "dftd3_version=4" in src
src=src.replace("  dftd3_version=4\n","  dftd3_version=4\n"+add,1)
# bump prefix
src=re.sub(r"prefix='[^']*'",lambda m:m.group().replace(m.group()[8:-1],m.group()[8:-1].replace("d3only","hse06d3")),src,1)
open(sys.argv[2],"w").write(src)
print("wrote",sys.argv[2])
