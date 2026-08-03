#!/usr/bin/env python3
"""Rebuild gd_neb/neb.in with a CONSISTENT molecular->chain path on the activated
Gd framework. Fixes the ordering bug that killed neb.x:
  - match the 8 I atoms between the molecular (-64.5) and chain (-132) endpoints
    by min-image assignment, reorder the chain endpoint to the molecular convention;
  - IDPP-interpolate 5 intermediates (collision-free, better than linear);
  - freeze the framework (reset to the exact activated positions) in every image;
  - reuse the existing neb.in's parse-validated &PATH/&CONTROL template, replacing
    only the BEGIN_POSITIONS block.
Also self-checks the max per-step per-atom displacement (< alat/2 needed)."""
import numpy as np, ase.io, os, shutil
from scipy.optimize import linear_sum_assignment
from ase.mep.neb import NEB

BASE = "/workspace/dimitry/calc/qe/i2_mof"
mol = ase.io.read(f"{BASE}/08_Gd/activated_i2/POSCAR_final_optimized")     # molecular, -64.5
cha = ase.io.read(f"{BASE}/08_Gd/activated_chain/POSCAR_final_optimized")  # chain,     -132
cell = mol.cell.array
sym = np.array(mol.get_chemical_symbols())
fwi = np.where(sym != "I")[0]; Ii = np.where(sym == "I")[0]
inv = np.linalg.inv(cell)
def mic(d):
    f = d @ inv; f -= np.round(f); return f @ cell
mI = mol.positions[Ii]; cI = cha.positions[Ii]
# assignment: mol_i -> chain_{perm[i]} minimizing sum of min-image dist^2
C = np.array([[(mic(mI[i]-cI[j])**2).sum() for j in range(8)] for i in range(8)])
_, perm = linear_sum_assignment(C)
cha2 = cha.copy()
cha2.positions[Ii] = cha.positions[Ii][perm]   # chain endpoint in molecular's I-order
# IDPP interpolate 5 intermediates
images = [mol] + [mol.copy() for _ in range(5)] + [cha2]
neb = NEB(images)
try:
    neb.interpolate(method="idpp"); method = "idpp"
except Exception as e:
    neb.interpolate(method="linear"); method = f"linear (idpp failed: {e!r})"
# freeze framework exactly in every image (set_positions, not in-place on a copy)
for im in images:
    p = im.get_positions()
    p[fwi] = mol.positions[fwi]
    im.set_positions(p)
fwmax = max(np.abs(im.positions[fwi] - mol.positions[fwi]).max() for im in images)
print(f"post-freeze framework max deviation from activated: {fwmax:.2e} A (want ~0)")
# self-check: max per-step per-atom min-image displacement
maxstep = 0.0; worst=None
for a in range(len(images)-1):
    d = images[a+1].positions - images[a].positions
    dw = mic(d); ms = np.linalg.norm(dw, axis=1)
    k = int(ms.argmax())
    if ms[k] > maxstep: maxstep=ms[k]; worst=(a,a+1,k)
print(f"interpolation: {method}")
print(f"max per-step per-atom disp: {maxstep:.2f} A  (alat/2=4.54 A; worst at img{worst[0]}->{worst[1]} atom#{worst[2]+1})")
# write neb.in (reuse template, replace positions)
nebf = f"{BASE}/gd_neb/neb.in"
if not os.path.exists(nebf+".bak"):
    shutil.copy(nebf, nebf+".bak")
tmpl = open(nebf+".bak").read()
pre = tmpl[:tmpl.index("BEGIN_POSITIONS")+len("BEGIN_POSITIONS")]
suf = tmpl[tmpl.index("END_POSITIONS"):]
out = [pre, ""]
labels = ["FIRST_IMAGE"] + ["INTERMEDIATE_IMAGE"]*5 + ["LAST_IMAGE"]
for lab, im in zip(labels, images):
    out.append(lab); out.append("ATOMIC_POSITIONS crystal")
    fr = im.get_scaled_positions(wrap=False)
    for s, f in zip(sym, fr):
        ifpos = "0 0 0" if s != "I" else "1 1 1"
        out.append(f"{s:<3s} {f[0]: .12f} {f[1]: .12f} {f[2]: .12f} {ifpos}")
out += ["", suf]
open(nebf, "w").write("\n".join(out))
print(f"wrote {nebf}  (old backed up to neb.in.bak)")

# ALSO regenerate pw_1..7.in (multi-file style) from the SAME images, so both
# input styles carry identical positions (no ambiguity which neb.x reads).
pw1 = f"{BASE}/gd_neb/pw_1.in"
if not os.path.exists(pw1+".bak"):
    shutil.copy(pw1, pw1+".bak")
tmpl = open(pw1+".bak").read()
hdr = tmpl[:tmpl.index("ATOMIC_POSITIONS")] + "ATOMIC_POSITIONS crystal\n"
for n, im in enumerate(images, start=1):
    lines = [hdr.rstrip("\n")]
    fr = im.get_scaled_positions(wrap=False)
    for s, f in zip(sym, fr):
        ifpos = "0 0 0" if s != "I" else "1 1 1"
        lines.append(f"{s:<4s} {f[0]:.12f} {f[1]:.12f} {f[2]:.12f}  {ifpos}")
    open(f"{BASE}/gd_neb/pw_{n}.in", "w").write("\n".join(lines) + "\n")
print(f"wrote pw_1..7.in (template backed up to pw_1.in.bak)")
