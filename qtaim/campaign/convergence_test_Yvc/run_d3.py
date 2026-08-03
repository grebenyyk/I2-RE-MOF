#!/usr/bin/env python3
"""Evaluate D3 dispersion via QE's bundled dftd3 (test.x) on a QE input geometry.
Usage: run_d3.py <qe_input.in> <func> <version>   -> prints EDISP_RY (Ry)."""
import sys, re, subprocess, os
import numpy as np
BOHR = 1.8897261254585  # Angstrom -> bohr
inp = sys.argv[1]; func = sys.argv[2]; ver = sys.argv[3]
TESTX = "/workspace/qe-7.5/dft-d3/test.x"
NVHPC = "/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11"
MPI = "/workspace/openmpi"
env = os.environ.copy()
env["PATH"] = MPI + "/bin:" + NVHPC + "/compilers/bin:" + env.get("PATH", "")
env["LD_LIBRARY_PATH"] = (MPI + "/lib:" + NVHPC + "/compilers/lib:" +
                          NVHPC + "/math_libs/lib:" + env.get("LD_LIBRARY_PATH", ""))
t = open(inp).read()
cellA = np.array([[float(x) for x in l.split()[:3]]
                  for l in t.split("CELL_PARAMETERS")[1].splitlines()[1:4]])
apos = []
for l in t.split("ATOMIC_POSITIONS")[1].split("K_POINTS")[0].splitlines():
    m = re.match(r"\s*([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", l)
    if m:
        apos.append((m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4))))
cellB = cellA * BOHR
lines = ["1.0000000000"]                       # a0 = 1 (latvecs already in bohr)
for i in range(3):
    lines.append("%.10f %.10f %.10f" % (cellB[i, 0], cellB[i, 1], cellB[i, 2]))
lines.append(str(len(apos)))
for sym, f1, f2, f3 in apos:
    cart = np.array([f1, f2, f3]) @ cellB
    lines.append("%-2s %.10f %.10f %.10f" % (sym, cart[0], cart[1], cart[2]))
stdin = "\n".join(lines) + "\n"
r = subprocess.run([TESTX, func, ver], input=stdin, capture_output=True, text=True, env=env)
m = re.search(r"EDISP_RY\s+(-?\d+\.\d+)", r.stdout + r.stderr)
print(float(m.group(1)) if m else "FAIL: " + (r.stdout[-400:] + r.stderr[-400:]))
