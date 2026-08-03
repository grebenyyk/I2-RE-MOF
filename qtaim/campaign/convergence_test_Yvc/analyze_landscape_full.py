#!/usr/bin/env python3
"""Full landscape: per variant, lowest molecular vs lowest polyiodide (HSE/4x4x2).
Scans all .out files for relevant tags. stabilization=(lowest_mol-lowest_poly)/4*1312.75."""
import re,os,glob
D="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/convergence_test_Yvc"
RYKJ=1312.75; N=4
def E(tag):
    f=f"{D}/{tag}.out"
    if not os.path.exists(f): return None
    m=re.findall(r"^!\s*total energy\s*=\s*(-?\d+\.\d+)",open(f,errors="ignore").read(),re.M)
    return float(m[-1]) if m else None
print("="*70+"\nFull landscape (HSE/4x4x2): lowest molecular vs lowest polyiodide\n"+"="*70)
for v in ["act","vc","st"]:
    print(f"\n=== Y-{v} ===")
    # polyiodide: authoritative + perturbed
    poly={}
    pe=E(f"rd_poly_Y{v}")
    if pe: poly["CIF (auth)"]=pe
    for t in sorted(glob.glob(f"{D}/*poly*Y{v}*.out"))+([f"{D}/hmol_Y{v}_S1.out"] if os.path.exists(f"{D}/hmol_Y{v}_S1.out") else []):
        tag=os.path.basename(t).replace(".out",""); e=E(tag)
        if e is not None and tag!=f"rd_poly_Y{v}": poly[tag]=e
    # molecular: per-variant S1 + S2/S3 (NOT cross-variant)
    mol={}
    s1={"act":"rd_mol_Yact","vc":"hmol_Yvc","st":None}[v]  # st S1 failed (SCF)
    if s1:
        e=E(s1)
        if e is not None: mol[s1]=e
    for pat in [f"hse_lmol_Y{v}_S*",f"hse_mmol_Y{v}"]:
        for t in sorted(glob.glob(f"{D}/{pat}.out")):
            tag=os.path.basename(t).replace(".out",""); e=E(tag)
            if e is not None: mol[tag]=e
    print(f"  Polyiodide configs ({len(poly)}):")
    for k,e in sorted(poly.items(),key=lambda x:x[1]):
        print(f"    {k:30s} {e:.5f}")
    print(f"  Molecular configs ({len(mol)}):")
    for k,e in sorted(mol.items(),key=lambda x:x[1]):
        print(f"    {k:30s} {e:.5f}")
    if poly and mol:
        lp=min(poly.values()); lm=min(mol.values())
        stab=(lm-lp)/N*RYKJ
        lp_name=min(poly,key=poly.get); lm_name=min(mol,key=mol.get)
        print(f"  => lowest poly: {lp_name} ({lp:.5f})")
        print(f"     lowest mol:  {lm_name} ({lm:.5f})")
        print(f"     stabilization = {stab:+.1f} kJ/I2 -> {'POLYIODIDE favored' if stab>0 else 'MOLECULAR favored'}")
