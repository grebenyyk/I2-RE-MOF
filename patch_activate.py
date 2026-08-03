#!/usr/bin/env python3
"""Patch the 7 activation inputs for robust strained-framework relaxation:
nosym=.true. (avoids checkallsym failures), mixing_beta 0.2->0.1 (conservative),
SCF conv_thr -> 1e-5 (strained frameworks stall at tighter), smaller BFGS trust radius."""
import os, glob
SYSTEMS = ["13_Y", "06_Tb", "12_Tb", "05_Cd", "07_Cd", "10_Cd", "17_Cd"]
for s in SYSTEMS:
    f = os.path.join(s, "activate_empty", "relax.in")
    t = open(f).read()
    # nosym
    if "nosym" not in t:
        t = t.replace("  nspin = 1", "  nspin = 1\n  nosym = .true.")
    # mixing 0.2 -> 0.1
    t = t.replace("mixing_beta = 0.2", "mixing_beta = 0.1")
    # SCF conv_thr -> 1e-5 (handle both 1e-6 and 1e-4 if present)
    for old in ["conv_thr = 1.0d-6", "conv_thr = 1.0d-4", "conv_thr = 1.d-6", "conv_thr = 1.d-4"]:
        t = t.replace(old, "conv_thr = 1.0d-5")
    # BFGS trust radius (prevent divergence)
    if "trust_radius" not in t:
        t = t.replace("  ion_dynamics = 'bfgs'",
                      "  ion_dynamics = 'bfgs'\n  trust_radius_max = 0.3\n  trust_radius_min = 1.0d-3\n  trust_radius_init = 0.1")
    open(f, "w").write(t)
    # clean stale output
    out = os.path.join(s, "activate_empty", "relax.out")
    if os.path.exists(out):
        os.remove(out)
    print(f"{s}: nosym={('nosym' in t)} mixing={('mixing_beta = 0.1' in t)} conv1e-5={('conv_thr = 1.0d-5' in t)} trust={('trust_radius' in t)}")
print("done")
