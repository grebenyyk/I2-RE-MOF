#!/usr/bin/env python3
"""Perturbed-seed test: are Y/Tb trapped above a deeper basin with short I-O contacts?

Take each structure's OWN optimized uniform_chain geometry, nudge the metal-coordinating I
atoms ~0.2 A toward their nearest framework O (seeding the 2.77 A I-O contact that Gd found),
then re-relax.  If the energy drops below the plain uniform_chain, Y/Tb were trapped; if they
relax back to their original I-O distance, the short contact is incompatible with their pore.

No Gd->Y mapping needed (uses each structure's own framework).
"""
import os, numpy as np, ase.io, ase
import build_inputs as bi
BASE = "/workspace/dimitry/calc/qe/i2_mof"
STRUCTS = [("13_Y", "Y"), ("06_Tb", "Tb")]
BUMP = 0.20     # A, displace coordinating I toward nearest O
MIND = 2.55     # A, don't push I closer than this to O (avoid severe clash)


def main():
    os.chdir(BASE)
    for s, metal in STRUCTS:
        a = ase.io.read(f"{s}/uniform_chain/POSCAR_final")   # optimized chain
        cell = np.array(a.cell); inv = np.linalg.inv(cell); pos = a.positions
        sy = np.array(a.get_chemical_symbols())
        def mic(i, j):
            d = pos[i]-pos[j]; d = d-cell.T@np.round(inv@d); return d, np.linalg.norm(d)
        Ii = [i for i, x in enumerate(sy) if x == "I"]
        Oi = [i for i, x in enumerate(sy) if x == "O"]
        Mi = [i for i, x in enumerate(sy) if x == metal]
        newpos = pos.copy()
        nudged = 0
        for i in Ii:
            # coordinating I = nearest metal < 3.5
            dM = min(np.linalg.norm(mic(i, j)[0]) for j in Mi)
            if dM > 3.6:
                continue
            # nearest framework O (min-image vector)
            best = min((mic(i, j) for j in Oi), key=lambda t: t[1])
            dvec, dd = best
            if dd < 1e-3:
                continue
            step = min(BUMP, max(0.0, dd - MIND))   # don't overshoot below MIND
            newpos[i] = pos[i] - step * dvec / dd
            nudged += 1
        bumped = ase.Atoms(symbols=sy, positions=newpos, cell=cell, pbc=True)
        outdir = f"{s}/uniform_chain_bump"
        os.makedirs(outdir, exist_ok=True)
        seed = f"{outdir}/POSCAR_seed"
        ase.io.write(seed, bumped, format="vasp", direct=True, vasp5=True)
        # clash check
        bp = bumped.positions; inv2 = np.linalg.inv(cell)
        def mic2(i, j):
            d = bp[i]-bp[j]; d = d-cell.T@np.round(inv2@d); return np.linalg.norm(d)
        fw = [i for i, x in enumerate(sy) if x != "I"]
        minIf = min(mic2(i, j) for i in Ii for j in fw)
        minII = min(mic2(i, j) for i in Ii for j in Ii if i > j)
        bi.build(seed, outdir, f"{s}_uniform_chain_bump")
        print(f"{s}: nudged {nudged} coordinating I toward O by {BUMP} A; min I-fw={minIf:.2f} min I-I={minII:.2f}")


if __name__ == "__main__":
    main()
