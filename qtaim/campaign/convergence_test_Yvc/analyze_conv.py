#!/usr/bin/env python3
"""Y-vc HSE convergence sweep analysis. For each knob, show energy vs setting and
the delta vs the finest (baseline = 4x4x2 / nqx 4 4 2 / ecutfock 150 / 50-400).
Stable to ~2 kJ/I2-equivalent (~0.0006 Ry/cell) = converged."""
import re,os
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
RYKJ=1312.75
def E(tag):
    f=f"{D}/{tag}.out"
    if not os.path.exists(f): return None
    m=re.findall(r"^!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f,errors="ignore").read(),re.M)
    return float(m[-1]) if m else None
base=E("conv_yvc_k4")
print(f"baseline (4x4x2 / nqx442 / ef150 / 50-400) = {base} Ry\n" if base else "baseline missing\n")
def row(tag,label):
    e=E(tag)
    if e is None: print(f"  {label:22s} (pending)"); return
    d=(e-base) if base else 0
    print(f"  {label:22s} E={e:.6f}  dE_vs_base={d*1000:+.2f} mRy  ({d*RYKJ/4:+.2f} kJ/I2)")
print("=== k-mesh scan (nqx follows k, ef150, 50/400) ===")
for t,l in [("conv_yvc_kG","k = 1x1x1 (Gamma)"),("conv_yvc_k2","k = 2x2x1"),("conv_yvc_k4","k = 4x4x2  <--base")]: row(t,l)
print("\n=== nqx scan (k=4x4x2 fixed, ef150, 50/400) ===")
for t,l in [("conv_yvc_nqx1","nqx = 1x1x1"),("conv_yvc_nqx2","nqx = 2x2x1"),("conv_yvc_k4","nqx = 4x4x2  <--base")]: row(t,l)
print("\n=== ecutfock scan (k4x4x2, nqx442, 50/400) ===")
for t,l in [("conv_yvc_ef100","ecutfock=100"),("conv_yvc_k4","ecutfock=150  <--base"),("conv_yvc_ef200","ecutfock=200")]: row(t,l)
print("\n=== ecutwfc/rho scan (k4x4x2, nqx442, ef150) ===")
for t,l in [("conv_yvc_ec40","ew/er=40/320"),("conv_yvc_k4","ew/er=50/400  <--base"),("conv_yvc_ec60","ew/er=60/480")]: row(t,l)
