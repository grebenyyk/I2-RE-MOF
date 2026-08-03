#!/usr/bin/env python3
"""Decomposed I2 sorption energies on ACTIVATED pore: dEtot = dEdisp + dEelec per I2.
Handles extension outputs (loaded_ext.out for the 4 nstep-capped chain jobs),
Gd (chain from activated_chain, mol from activated_i2), and the rest.
Usage: python3 collect_activated.py  -> prints CHAIN + MOLECULAR tables."""
import re, os
EI2=-757.79989474; ED3_I2=-0.00354048; RY=13.605693009
BASE="/workspace/dimitry/calc/qe/i2_mof"
def last(fn,pat):
    if not os.path.exists(fn): return None
    m=re.findall(pat,open(fn,errors="ignore").read()); return float(m[-1]) if m else None
def Et(fn): return last(fn,r"!\s*total energy\s*=\s*(-?\d+\.\d+)")
def Ed(fn): return last(fn,r"DFT-D3 Dispersion\s*=\s*(-?\d+\.\d+)")
SYS=[("08_Gd",4),("06_Tb",4),("13_Y",4),("12_Tb",6),("17_Cd",2),("07_Cd",6),("10_Cd",5),("05_Cd",5)]
def chain_src(s):
    if s=="08_Gd": return f"{BASE}/08_Gd/activated_chain/loaded.out"
    ext=f"{BASE}/{s}/activated_i2/loaded_ext.out"
    if Et(ext) is not None: return ext
    return f"{BASE}/{s}/activated_i2/loaded.out"
def mol_src(s):
    if s=="08_Gd": return f"{BASE}/08_Gd/activated_i2/loaded.out"  # Gd molecular = -64.5
    return f"{BASE}/{s}/activated_i2_mol/loaded.out"
def empty_src(s): return f"{BASE}/{s}/activate_empty/relax.out"
def table(srcfn,label):
    rows=[]
    for s,n in SYS:
        eL=Et(srcfn(s)); dL=Ed(srcfn(s)); eE=Et(empty_src(s)); dE=Ed(empty_src(s))
        if eL and eE and dL and dE:
            et=(eL-eE-n*EI2)*RY/n; ed=(dL-dE-n*ED3_I2)*RY/n
            rows.append((s,n,et,et*96.485,ed*96.485,(et-ed)*96.485))
        else: rows.append((s,n,None,None,None,None))
    got=[r for r in rows if r[2] is not None]
    if not got: print(f"\n=== {label}: none complete yet ===\n"); return
    print(f"\n=== {label} on activated pore ===")
    print(f"{'System':<8}{'NI2':>4}{'dEtot eV':>10}{'dEtot kJ':>10}{'dEdisp kJ':>11}{'dEelec kJ':>11}")
    for s,n,et,ev,ed,ee in rows:
        if et is None: print(f"{s:<8}{n:>4}  (pending)"); continue
        print(f"{s:<8}{n:>4}{et:>10.4f}{ev:>10.2f}{ed:>11.2f}{ee:>11.2f}")
table(chain_src,"CHAIN")
table(mol_src,"MOLECULAR")
