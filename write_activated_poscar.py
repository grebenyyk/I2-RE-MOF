#!/usr/bin/env python3
"""Write VASP POSCARs of the activated (relaxed-empty) frameworks for every
system whose relax.out has reached JOB DONE. Reads the FINAL QE image via ASE
(cell in A, final fractional coords), groups atoms by species in first-appearance
order. Idempotent: re-run anytime; only completed systems are written."""
import ase.io, os, re
BASE = "/workspace/dimitry/calc/qe/i2_mof"
SYSTEMS = ["13_Y", "06_Tb", "12_Tb", "05_Cd", "07_Cd", "10_Cd", "17_Cd"]
for s in SYSTEMS:
    out = f"{BASE}/{s}/activate_empty/relax.out"
    if not os.path.isfile(out) or "JOB DONE" not in open(out, errors="ignore").read():
        print(f"{s}: not done yet, skip")
        continue
    a = ase.io.read(out, format="espresso-out", index=-1)
    cell = a.cell.array
    frac = a.get_scaled_positions(wrap=False)
    syms = a.get_chemical_symbols()
    order, cnt = [], []
    for sp in syms:
        if sp not in order:
            order.append(sp); cnt.append(1)
        else:
            cnt[order.index(sp)] += 1
    L = [f"POSCAR_activated_empty  ({s} relaxed-empty framework, from relax.out final coords)",
         "1.0"]
    for r in cell:
        L.append(f"   {r[0]: .10f}   {r[1]: .10f}   {r[2]: .10f}")
    L.append("   " + "   ".join(f"{sp:>4s}" for sp in order))
    L.append("   " + "   ".join(f"{c:4d}" for c in cnt))
    L.append("Direct")
    for sp in order:
        for s2, f in zip(syms, frac):
            if s2 == sp:
                L.append(f"   {f[0]: .12f}   {f[1]: .12f}   {f[2]: .12f}")
    dst = f"{BASE}/{s}/activate_empty/POSCAR_activated_empty"
    open(dst, "w").write("\n".join(L) + "\n")
    print(f"{s}: {len(syms)} atoms  species={dict(zip(order,cnt))}  -> {dst}")
