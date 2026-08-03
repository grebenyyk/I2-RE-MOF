#!/usr/bin/env python3
"""HSE Born-Haber decomposition (corrected cycle), 3 Y variants.
deltaE_chain-assembly = (E*_I - 4*E_gas)/4
deltaE_framework-deform = (E*_host - E0_host)/4   [0 for st/act; vc uses dc_hostld_Yvc]
deltaE_anchoring = (E_loaded - E*_host - E*_I)/4
sum = deltaE_adsorption (checks vs Table 1)."""
import re,os
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
RYKJ=1312.75; N=4
def E(tag):
    f=f"{D}/{tag}.out"
    if not os.path.exists(f): return None
    m=re.findall(r"^!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f,errors="ignore").read(),re.M)
    return float(m[-1]) if m else None
gas=E("rd_gas")
print("="*60+"\nHSE Born-Haber decomposition (Y, kJ/mol per I2)\n"+"="*60)
print(f"{'System':12}{'chain-assm':>11}{'deform':>9}{'anchor':>9}{'sum':>9}{'Table1':>9}")
print("-"*58)
for v in ["st","act","vc"]:
    L=E(f"rd_poly_Y{v}"); E0=E(f"rd_empty_Y{v}"); I=E(f"dc_isubl_Y{v}")
    Hstar = E("dc_hostld_Yvc") if v=="vc" else E0   # st/act: Host*=Host0
    if None in (L,E0,I,gas) or (v=="vc" and Hstar is None):
        print(f"Y-{v:3}  (missing)"); continue
    chain=(I-N*gas)/N*RYKJ
    deform=(Hstar-E0)/N*RYKJ
    anchor=(L-Hstar-I)/N*RYKJ
    s=chain+deform+anchor
    t1=(L-E0-N*gas)/N*RYKJ
    print(f"I2@Y-MOF{v:<4}{chain:>11.2f}{deform:>9.2f}{anchor:>9.2f}{s:>9.2f}{t1:>9.2f}")
