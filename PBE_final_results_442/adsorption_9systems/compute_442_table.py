#!/usr/bin/env python3
"""Compute the 4x4x2-relaxed adsorption-energy table for the 9 systems.
Reads rd_poly_442_<sys>.out (loaded), rd_empty_<sys>.out (framework), rd_gas.out (gas I2),
all at degauss=0.02, 4x4x2 (loaded) / Gamma (empty, gas). dEads per I2 = (E_loaded -
E_empty - 4*E_gas)/4. Only uses energies from CONVERGED outputs; marks others PENDING.
Run after SP02 (clean 0.02 SPs on the 0.04-relaxed act/st systems) completes."""
import re,glob,os,shutil
D="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results_442/adsorption_9systems"
G="/workspace/dimitry/calc/qe/i2_mof/PBE_final_results/adsorption_9systems"  # Γ-relaxed (Table 1) for comparison
REL="/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/relax_442"
RY2KJ=13.6057*96.485
SYS=["Yvc","Gdvc","Tbvc","Yact","Gdact","Tbact","Yst","Gdst","Tbst"]
# --- refresh: pull the latest clean 0.02 SP outputs (SP02) for the systems relaxed at 0.04 ---
# vc + Yact were relaxed at 0.02 -> their rd_poly_442 .out is the relaxation output (already clean).
# Gdact/Tbact/Yst/Gdst/Tbst were relaxed at 0.04 -> use the clean 0.02 SP02 single-points.
SP02SRC={"Gdact":f"{REL}/sp02/sp02_Gdact.out","Tbact":f"{REL}/sp02/sp02_Tbact.out",
         "Yst":f"{REL}/sp02/sp02_Yst.out","Gdst":f"{REL}/sp02/sp02_Gdst.out","Tbst":f"{REL}/sp02/sp02_Tbst.out"}
def is_conv(f):
    if not os.path.exists(f): return False
    s=open(f,errors="ignore").read()
    return "convergence has been achieved" in s and "!" in s
for s,src in SP02SRC.items():
    insrc=src.replace(".out",".in")
    if os.path.exists(insrc):
        shutil.copy(insrc,f"{D}/rd_poly_442_{s}.in")  # real SP02 input
    if is_conv(src):
        shutil.copy(src,f"{D}/rd_poly_442_{s}.out")  # converged SP02 output
def Econ(f):
    if not os.path.exists(f): return None,None
    s=open(f,errors="ignore").read()
    conv = ("bfgs converged" in s or "Damped Dynamics: convergence achieved" in s
            or "convergence has been achieved" in s)
    m=re.findall(r"!\s*total energy\s*=\s*(-?\d+\.\d+)",s)
    return (float(m[-1]) if m else None, conv)
gas,_=Econ(f"{D}/rd_gas.out")
empty={s:Econ(f"{D}/rd_empty_{s}.out")[0] for s in SYS}
rows=[]
for s in SYS:
    el,c=Econ(f"{D}/rd_poly_442_{s}.out")
    eg,cg=Econ(f"{G}/rd_poly_{s}.out")  # Γ-relaxed production energy (Table 1)
    if el is not None and c and gas is not None and empty[s] is not None:
        da442=(el-empty[s]-4*gas)*RY2KJ/4.0
        dag=(eg-empty[s]-4*gas)*RY2KJ/4.0 if eg is not None else None
        rows.append((s,da442,dag,da442-dag if dag is not None else None))
    else:
        rows.append((s,None,None,None))
# write csv
with open(f"{D}/adsorption_table_442.csv","w") as fp:
    fp.write("system,dEads_442(kJ/I2),dEads_Gamma(kJ/I2),delta(442-Gamma)\n")
    for s,a,b,d in rows:
        fp.write(f"{s},{a if a is not None else 'PENDING'},{b if b is not None else ''},{d if d is not None else ''}\n")
print(f"{'system':7} {'dEads 4x4x2-relax':>18} {'dEads Gamma(prod)':>18} {'delta':>9}")
print("-"*56)
for s,a,b,d in rows:
    print(f"{s:7} {(f'{a:.2f}' if a is not None else 'PENDING'):>18} {(f'{b:.2f}' if b is not None else '-'):>18} {(f'{d:+.2f}' if d is not None else '-'):>9}")
print(f"\n(gas={gas:.4f}, RY2KJ={RY2KJ:.3f}; empty/gas at Gamma 0.02; loaded at 4x4x2 0.02)")
