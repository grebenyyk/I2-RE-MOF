#!/usr/bin/env python3
"""Build 24x24 I-sublattice Hessian -> frequencies, for a given raman dir.
Usage: python3 hessian_442.py <dir>  (dir has tags.txt + disp_<tag>.out)
Same logic as raman_driver.build_hessian: Phi[:,col]=-(F+ - F-)/(2*delta_bohr),
mass-weight m_I, RY_TO_CM=3634.77."""
import sys, re, numpy as np
M_I=126.904; RY_TO_CM=3634.77; DELTA_A=0.01; N_I=8; bohr=0.529177210563
d=sys.argv[1] if len(sys.argv)>1 else "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/Y_deep/kmesh_test/raman_442"
def parse_forces(tag):
    o=f"{d}/disp_{tag}.out"; t=open(o,errors="ignore").read()
    # forces on atoms 73..80 (the 8 I atoms), Ry/bohr
    rows=[]
    for atm in range(73,81):
        m=re.search(rf"atom\s+{atm}\s+type\s+\d+\s+force\s*=\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",t)
        if not m: return None
        rows.append([float(m.group(i)) for i in (1,2,3)])
    return np.array(rows)  # (8,3)
raw=[l.split() for l in open(f"{d}/tags.txt") if l.strip()]
tags=[(int(r[0]),int(r[1]),r[2],r[3]) for r in raw]
Phi=np.zeros((3*N_I,3*N_I)); delta_bohr=DELTA_A/bohr; F={}
miss=[]
for i,a,sign,tag in tags:
    f=parse_forces(tag)
    if f is None: miss.append(tag); continue
    F[(i,a,sign)]=f.flatten()
for i,a,_,_ in tags:
    col=3*i+a
    if (i,a,'p') in F and (i,a,'n') in F:
        Phi[:,col]=-(F[(i,a,'p')]-F[(i,a,'n')])/(2*delta_bohr)
    else: miss.append(f"col{col}")
D=Phi/M_I  # mass-weight: sqrt(m_i*m_j)=M_I (all 8 atoms are I)
evals,evecs=np.linalg.eigh(D)
freqs=np.sign(evals)*np.sqrt(np.abs(evals))*RY_TO_CM
order=np.argsort(freqs)[::-1]
label=d.rstrip('/').split('/')[-1]
print(f"=== {label}: 24 I-sublattice frequencies (cm-1) ===")
print(f"{'mode':>4} {'freq':>9}")
for k,ix in enumerate(order):
    flag=" (imag)" if freqs[ix]<0 else ""
    print(f"{k+1:>4} {freqs[ix]:>9.2f}{flag}")
print(f"real>0: {(freqs>0).sum()} | imag: {(freqs<0).sum()} | missing: {len(set(miss))}")
np.save(f"{d}/freqs_442.npy", freqs[order]) if "442" in label else None
