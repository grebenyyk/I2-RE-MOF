#!/usr/bin/env python3
"""Verify-suite analysis: (1) vc k-shift at 60/800 (clean, vs the 60/800 Gamma refs),
(2) MP-vs-mv smearing check, (3) molecular-vs-polyiodide at 4x4x2."""
import re,os
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
BASE="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"; ROOT="/workspace/dimitry/calc/qe/i2_mof"
RY=1312.75
def E(f):
    if not os.path.exists(f): return None
    m=re.findall(r"^!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f,errors="ignore").read(),re.M)
    return float(m[-1]) if m else None
GAS=E(f"{ROOT}/i2_vac_ref/scf.out")
ACT_EMPTY_Y=E(f"{ROOT}/13_Y/activate_empty/relax.out")  # Y-act empty (Gamma)
print("="*64+"\n(1) VC k-SHIFT at 60/800  (E(442)-E(Gamma))/4, kJ/I2  [clean: same cutoff]\n"+"="*64)
SCFL={"Y_deep":-5767.52921047,"Gd_deep":-6865.27509945,"Tb_deep":-6925.32970510}  # scf_loaded.out (60/800 Gamma)
SCFI={"Y_deep":-3031.42652513,"Gd_deep":-3031.42557816,"Tb_deep":-3031.42904301}  # scf_i2only.out (60/800 Gamma)
PROD_VC={"Y_deep":-163.10,"Gd_deep":-167.38,"Tb_deep":-163.20}
print(f"{'system':9}{'loaded shift':>13}{'isubl shift':>13}{'vc dEads corr':>14}")
for v in ["Y_deep","Gd_deep","Tb_deep"]:
    L442=E(f"{D}/verify_vcL_{v}.out"); I442=E(f"{D}/verify_vcI_{v}.out")
    if L442 is None: print(f"{v:9} (missing)"); continue
    sl=(L442-SCFL[v])/4*RY; si=(I442-SCFI[v])/4*RY if I442 else float('nan')
    print(f"{v:9}{sl:>13.1f}{si:>13.1f}{PROD_VC[v]+sl:>14.1f}")
print("\n"+"="*64+"\n(2) MP vs mv SMEARING check  (kJ/I2, should agree to <~3)\n"+"="*64)
# Y-act: redo442_loaded_Y_actc (mv, 50/400 442) vs verify_mp_yact (mp, 50/400 442)
mv_act=E(f"{D}/redo442_loaded_Y_actc.out"); mp_act=E(f"{D}/verify_mp_yact.out")
if mv_act and mp_act: print(f"Y-act (50/400 442): mv={mv_act:.5f}  mp={mp_act:.5f}  diff={(mp_act-mv_act)/4*RY:+.2f} kJ/I2  -> {'AGREE' if abs((mp_act-mv_act)/4*RY)<3 else 'DISAGREE'}")
# Y-vc: verify_vcL_Y_deep (mv, 60/800 442) vs verify_mp_yvc (mp, 60/800 442)
mv_vc=E(f"{D}/verify_vcL_Y_deep.out"); mp_vc=E(f"{D}/verify_mp_yvc.out")
if mv_vc and mp_vc: print(f"Y-vc  (60/800 442): mv={mv_vc:.5f}  mp={mp_vc:.5f}  diff={(mp_vc-mv_vc)/4*RY:+.2f} kJ/I2  -> {'AGREE' if abs((mp_vc-mv_vc)/4*RY)<3 else 'DISAGREE'}")
print("\n"+"="*64+"\n(3) MOLECULAR vs POLYIODIDE at 4x4x2 (Y-act, 50/400)  kJ/I2\n"+"="*64)
poly=E(f"{D}/redo442_loaded_Y_actc.out"); mol=E(f"{D}/verify_mol_yact.out")
if poly and mol and ACT_EMPTY_Y:
    dp=(poly-ACT_EMPTY_Y-4*GAS)/4*RY; dm=(mol-ACT_EMPTY_Y-4*GAS)/4*RY
    print(f"polyiodide dEads(442) = {dp:.1f}")
    print(f"molecular  dEads(442) = {dm:.1f}")
    print(f"polyiodide favored by {dp-dm:.1f} kJ/I2  -> {'POLYIODIDE stable' if dp<dm else 'MOLECULAR favored (!)'}")
    print(f"(at Gamma, polyiodide was favored by ~58 kJ/I2)")
