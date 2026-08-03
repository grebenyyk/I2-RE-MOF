#!/usr/bin/env python3
"""Definitive I2@MOF sorption table. ONE reference set for activated + strained
so the two columns are directly comparable. Strained dE uses strained-empty
(empty_mof.out) as baseline; activated dE uses activated-empty (activate_empty/
relax.out). Decomposition dEtot = dEdisp + dEelec, per I2."""
import re, os
BASE="/workspace/dimitry/calc/qe/i2_mof"
EI2=-757.79989474; ED3_I2=-0.00354048; RY=13.605693009; KJ=96.48533212
SYS=[("08_Gd",4),("06_Tb",4),("13_Y",4),("12_Tb",6),("17_Cd",2),("07_Cd",6),("10_Cd",5),("05_Cd",5)]
NAME={"08_Gd":"Gd-MOF","06_Tb":"Tb-MOF","13_Y":"Y-MOF","12_Tb":"12-Tb","17_Cd":"17-Cd",
     "07_Cd":"07-Cd","10_Cd":"10-Cd","05_Cd":"05-Cd"}
def last(fn,pat):
    if not os.path.exists(fn): return None
    m=re.findall(pat,open(fn,errors="ignore").read()); return float(m[-1]) if m else None
def Et(fn): return last(fn,r"!\s*total energy\s*=\s*(-?\d+\.\d+)")
def Ed(fn): return last(fn,r"DFT-D3 Dispersion\s*=\s*(-?\d+\.\d+)")
def calc(loaded, empty, n):
    eL,dL=Et(loaded),Ed(loaded); eE,dE=Et(empty),Ed(empty)
    et=(eL-eE-n*EI2)*RY/n*KJ; ed=(dL-dE-n*ED3_I2)*RY*KJ/n; ee=et-ed
    return et,ed,ee
# --- sources ---
ACT_CHAIN={"08_Gd":"08_Gd/activated_chain/loaded.out"}
for s in ["06_Tb","13_Y","12_Tb","17_Cd","05_Cd"]:
    p=f"{s}/activated_i2/loaded_ext.out"; ACT_CHAIN[s]=p if Et(f"{BASE}/{p}") else f"{s}/activated_i2/loaded.out"
for s in ["07_Cd","10_Cd"]: ACT_CHAIN[s]=f"{s}/activated_i2/loaded.out"
ACT_MOL={"08_Gd":"08_Gd/activated_i2/loaded.out"}
for s in ["06_Tb","13_Y","12_Tb","17_Cd","07_Cd","10_Cd","05_Cd"]: ACT_MOL[s]=f"{s}/activated_i2_mol/loaded.out"
ACT_EMPTY=lambda s:f"{s}/activate_empty/relax.out"
STR_CHAIN={"08_Gd":"08_Gd/uniform_chain/loaded.out","06_Tb":"06_Tb/uniform_chain/loaded.out","13_Y":"13_Y/uniform_chain/loaded.out",
           "12_Tb":"12_Tb/loaded_tight.out","17_Cd":"17_Cd/loaded.out","07_Cd":"07_Cd/loaded.out","10_Cd":"10_Cd/loaded.out","05_Cd":"05_Cd/loaded.out"}
STR_EMPTY=lambda s:f"{s}/empty_mof.out"
# --- print ---
def fmt(a,b): return f"{a:8.2f} ({b:8.2f})"  # activated (strained)
print("I2@MOF sorption energies - CHAIN motif")
print("format: activated value (strained value in brackets). Per I2. Single I2 ref (-757.79989474 Ry).")
print()
print(f"{'System':<8}{'NI2':>4}{'dEtot/I2 kJ/mol':>22}{'dEdisp/I2 kJ/mol':>22}{'dEelec/I2 kJ/mol':>22}")
print("-"*78)
for s,n in SYS:
    aet,aed,aee=calc(f"{BASE}/{ACT_CHAIN[s]}",f"{BASE}/{ACT_EMPTY(s)}",n)
    set_,seed,see=calc(f"{BASE}/{STR_CHAIN[s]}",f"{BASE}/{STR_EMPTY(s)}",n)
    print(f"{NAME[s]:<8}{n:>4}{fmt(aet,set_):>22}{fmt(aed,seed):>22}{fmt(aee,see):>22}")
print()
print("MOLECULAR motif (activated only):")
print(f"{'System':<8}{'NI2':>4}{'dEtot/I2 kJ/mol':>16}{'dEdisp/I2 kJ/mol':>17}{'dEelec/I2 kJ/mol':>17}")
for s,n in SYS:
    aet,aed,aee=calc(f"{BASE}/{ACT_MOL[s]}",f"{BASE}/{ACT_EMPTY(s)}",n)
    print(f"{NAME[s]:<8}{n:>4}{aet:>16.2f}{aed:>17.2f}{aee:>17.2f}")
print()
print(f"I2_vac ref: E={EI2} Ry, D3={ED3_I2} Ry (i2_vac_ref/scf.out). Strained = guest-imposed CoRE-MOF pore; empty baseline = empty_mof.out. Activated empty baseline = activate_empty/relax.out.")
