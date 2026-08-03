#!/usr/bin/env python3
"""Molecular landscape: for each variant, all molecular HSE energies vs polyiodide.
Lowest molecular = best molecular can do. stabilization=(E_mol-E_poly)/4*1312.75;
negative = molecular favored. act includes rd_mol_Yact (the -16.7 result)."""
import re,os
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
RYKJ=1312.75; N=4
def E(tag):
    f=f"{D}/{tag}.out"
    if not os.path.exists(f): return None
    m=re.findall(r"^!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f,errors="ignore").read(),re.M)
    return float(m[-1]) if m else None
print("="*62+"\nMolecular landscape vs polyiodide (HSE/4x4x2, kJ/I2)\n"+"="*62)
for v in ["st","act","vc"]:
    ep=E(f"rd_poly_Y{v}")
    mols=[]
    for s in range(1,6):
        e=E(f"hmol_Y{v}_S{s}")
        if e is not None: mols.append((f"seed{s}",e))
    if v=="act":
        e=E("rd_mol_Yact")
        if e is not None: mols.append(("rd_mol_Yact",e))
    if ep is None or not mols:
        print(f"\nY-{v}: (poly={'y' if ep else 'n'}, {len(mols)} molecular) pending"); continue
    print(f"\nY-{v}  (polyiodide E={ep:.5f}, {len(mols)} molecular configs):")
    best=min(mols,key=lambda x:x[1])
    for nm,e in sorted(mols,key=lambda x:x[1]):
        st=(e-ep)/N*RYKJ
        flag="  <-- lowest" if nm==best[0] else ""
        print(f"  {nm:14s} E={e:.5f}  stabilization={st:+7.1f}{flag}")
    bst=(best[1]-ep)/N*RYKJ
    print(f"  => lowest molecular {best[0]}: stabilization={bst:+.1f} -> "
          f"{'POLYIODIDE favored' if bst>0 else 'MOLECULAR favored'}")
