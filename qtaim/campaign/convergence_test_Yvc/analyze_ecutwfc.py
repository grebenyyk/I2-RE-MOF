#!/usr/bin/env python3
"""ecutwfc BINDING convergence (Y-vc) + 6x6x3 k-check. All geometries CIF-correct.
binding = (loaded - empty - 4*gas)/4 * 1312.75. loaded: conv_yvc_ec40/k4/ec60
(CIF). empty: chk_empty_yvc_ec{40,50,60} (CIF). gas: chk_gas_ec{40,60}+hse_i2_vac."""
import re,os
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
RYKJ=1312.75; N=4
def E(tag):
    f=f"{D}/{tag}.out"
    if not os.path.exists(f): return None
    m=re.findall(r"^!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f,errors="ignore").read(),re.M)
    return float(m[-1]) if m else None
LOADED={40:"conv_yvc_ec40",50:"conv_yvc_k4",60:"conv_yvc_ec60"}
EMPTY={40:"chk_empty_yvc_ec40",50:"chk_empty_yvc_ec50",60:"chk_empty_yvc_ec60"}
GAS={40:"chk_gas_ec40",50:"hse_i2_vac",60:"chk_gas_ec60"}
print("=== ecutwfc BINDING convergence (Y-vc) ===")
prev=None
for ew in [40,50,60]:
    L=E(LOADED[ew]); F=E(EMPTY[ew]); g=E(GAS[ew])
    miss=[k for k,v in [('L',L),('F',F),('g',g)] if v is None]
    if miss: print(f"  ew/er={ew}/{8*ew}: missing {miss}"); continue
    b=(L-F-N*g)/N*RYKJ
    d=f"  (delta_vs_prev={b-prev:+.2f} kJ/I2)" if prev is not None else ""
    print(f"  ew/er={ew}/{8*ew}: binding={b:.2f} kJ/I2{d}   [loaded={L:.5f} empty={F:.5f} gas={g:.5f}]")
    prev=b
print("\n=== 6x6x3 k-check (loaded, nqx=k, ef100, 50/400) ===")
k4=E("conv_yvc_k4"); k6=E("chk_yvc_k663")
if k6 is None: print("  6x6x3 pending")
else:
    d=(k6-k4)*RYKJ/N
    print(f"  4x4x2 = {k4:.6f}   6x6x3 = {k6:.6f}   delta={d:+.2f} kJ/I2  -> {'4x4x2 CONVERGED (<2 kJ)' if abs(d)<2 else '6x6x3 NEEDED (>'+f'{abs(d):.0f} kJ)'}")
