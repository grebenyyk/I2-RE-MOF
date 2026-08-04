#!/usr/bin/env python3
"""4x4x2-relaxed 9-system adsorption table with D3/PBE decomposition.
ΔEads = ΔEdisp(D3) + ΔEelec(PBE), all per I2 (4 I2/cell).
D3 energy prints directly as 'DFT-D3 Dispersion = X Ry' in every output.
Self-refreshing: copies converged SP02 (.in+.out) for the act/st systems relaxed at 0.04."""
import re,os,shutil
D="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results_442/adsorption_9systems"
G="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/adsorption_9systems"  # Γ-relaxed for comparison
REL="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/relax_442"
RY2KJ=13.6057*96.485; N=4
# refresh: pull real SP02 (.in+.out) for act/st systems relaxed at degauss=0.04
SP02SRC={"Gdact":f"{REL}/sp02/sp02_Gdact.out","Tbact":f"{REL}/sp02/sp02_Tbact.out",
         "Yst":f"{REL}/sp02/sp02_Yst.out","Gdst":f"{REL}/sp02/sp02_Gdst.out","Tbst":f"{REL}/sp02/sp02_Tbst.out"}
def is_conv(f):
    if not os.path.exists(f): return False
    s=open(f,errors="ignore").read()
    return "convergence has been achieved" in s and "!" in s
for s,src in SP02SRC.items():
    insrc=src.replace(".out",".in")
    if os.path.exists(insrc): shutil.copy(insrc,f"{D}/single_points/rd_poly_442_{s}.in")
    if is_conv(src): shutil.copy(src,f"{D}/single_points/rd_poly_442_{s}.out")
def parse(fn):
    """return (E_total, D3, converged) or (None,None,False)."""
    if not os.path.exists(fn): return None,None,False
    s=open(fn,errors="ignore").read()
    conv=("bfgs converged" in s or "Damped Dynamics: convergence achieved" in s or "convergence has been achieved" in s)
    et=re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)",s)
    d3=re.findall(r"DFT-D3 Dispersion\s*=\s*(-?\d+\.\d+)",s)
    return (float(et[-1]) if et else None, float(d3[-1]) if d3 else None, conv)
gasEt,gasD3,_=parse(f"{D}/rd_gas.out")
SYS=["Gdvc","Tbvc","Yvc","Gdact","Tbact","Yact","Gdst","Tbst","Yst"]  # printed grouped below
def decomp(s):
    el,dl,cl=parse(f"{D}/single_points/rd_poly_442_{s}.out")
    ee,de,_=parse(f"{D}/rd_empty_{s}.out")
    if None in (el,dl,ee,de,gasEt,gasD3) or not cl:
        return None
    dEads=(el-ee-N*gasEt)*RY2KJ/N
    dEdisp=(dl-de-N*gasD3)*RY2KJ/N
    return (dEads,dEdisp,dEads-dEdisp)
# print grouped (st, act, vc) x (Gd, Tb, Y) + write csv
rows=[]
print(f"{'System':14} {'ΔEads/I2':>10} {'ΔED3/I2':>10} {'ΔEPBE/I2':>10}   (kJ/mol)")
print("-"*52)
for var,grp in [("st","Solvent-templated"),("act","Activated"),("vc","Variable-cell")]:
    print(grp)
    for m in ["Gd","Tb","Y"]:
        s=f"{m}{var}"; r=decomp(s)
        if r: print(f"I2@{m}-MOF{var:<5} {r[0]:10.2f} {r[1]:10.2f} {r[2]:10.2f}"); rows.append((f"I2@{m}-MOF{var}",r[0],r[1],r[2]))
        else: print(f"I2@{m}-MOF{var:<5} {'PENDING':>10}"); rows.append((f"I2@{m}-MOF{var}",None,None,None))
import csv
with open(f"{D}/adsorption_table_442.csv","w",newline="") as fp:
    w=csv.writer(fp); w.writerow(["system","dEads_I2","dEdisp_I2","delec_I2","kJ/mol; delec=dEads-dEdisp"])
    for name,a,d,e in rows:
        w.writerow([name, f"{a:.2f}" if a is not None else "PENDING",
                    f"{d:.2f}" if d is not None else "", f"{e:.2f}" if e is not None else ""])
print(f"\n(gas Et={gasEt:.4f} D3={gasD3:.5f}; empty/gas at Gamma 0.02; loaded at 4x4x2 0.02; {N} I2/cell)")
