#!/usr/bin/env python3
"""CShM analysis of the 9 i2_mof structures (3 metals Y/Gd/Tb x {pristine, stripped, POSCAR}).

Protocol:
  pristine  : 1 unique metal (spglib P2_1/c), CN=8  -> 8-vertex references
  stripped  : 1 unique metal (spglib recovers P2_1/c from the lying P1 CIF), CN=6
  POSCAR    : true P1 -> all 4 metals; each scored at CN=6 (6 O, framework)
              AND CN=7 (6 O + inner iodine <=3.7 A)
Cutoffs: M-O/M-N = 2.75 A; M-I = 3.7 A (excludes 2nd-shell I ~4.0-4.3 A).
Donor-atom convention vs N-vertex references, central_atom=0.
"""
import warnings, numpy as np
warnings.filterwarnings("ignore")
from ase.io import read
from ase.neighborlist import neighbor_list
import spglib
from cosymlib import Geometry
from cosymlib.shape import Shape
from cosymlib.shape.tools import get_structure_references

DIR = "/workspace/dimitry/calc/qe/i2_mof/cshm"
CUT_O, CUT_I = 2.75, 3.7
DONOR_O = {"O", "N", "S", "P", "F", "Cl", "Br"}
METALS = set("Y Gd Tb La Ce Pr Nd Sm Eu Dy Ho Er Tm Yb Lu Sc Zr Hf Ti Zn".split())

LABEL_NAMES = {
    'OC-6': 'octahedron', 'TPR-6': 'trigonal prism', 'HP-6': 'hexagon',
    'PPY-6': 'pentagonal pyramid', 'JPPY-6': 'J pentagonal pyramid',
    'PBPY-7': 'pentagonal bipyramid', 'COC-7': 'capped octahedron',
    'CTPR-7': 'capped trig prism', 'HPY-7': 'heptagonal pyramid', 'HP-7': 'heptagon',
    'JPBPY-7': 'J pentagonal bipyramid', 'JETPY-7': 'J edge trig pyramid',
    'SAPR-8': 'square antiprism', 'CU-8': 'cube', 'TDD-8': 'triang dodecahedron',
    'BTPR-8': 'bicapped trig prism', 'HBPY-8': 'hexagonal bipyramid',
    'OP-8': 'octagon prism', 'HPY-8': 'heptagonal pyramid',
    'JGBF-8': 'J gyro bifastigium', 'JETBPY-8': 'J edge bicapped trig prism',
    'JBTPR-8': 'J bicapped trig prism', 'JSD-8': 'J snub disphenoid',
    'TT-8': 'truncated tetrahedron', 'ETBPY-8': 'edge trig bicapped pyramid',
}
def desc(l): return LABEL_NAMES.get(l, l)

def donors_of(at, idx, include_I):
    """(dist, j, sym, vec) for coordinating donors of atom idx across PBC.
    vec is the metal->donor minimum-image vector, so metal_pos+vec is the
    donor's true position relative to the metal (correct across cell edges)."""
    i_a, j_a, d_a, v_a = neighbor_list("ijdD", at, CUT_I + 0.05)
    out = []
    for d, j, v in zip(d_a[i_a == idx], j_a[i_a == idx], v_a[i_a == idx]):
        if j == idx:
            continue
        s = at[j].symbol
        if s in DONOR_O and d <= CUT_O:
            out.append((float(d), int(j), s, np.asarray(v)))
        elif s == "I" and include_I and d <= CUT_I:
            out.append((float(d), int(j), s, np.asarray(v)))
    out.sort(key=lambda x: x[0])
    return out

def cshm_rows(donor_abs, cn):
    refs = [str(r) for r in get_structure_references(cn)]
    coords = donor_abs - donor_abs.mean(axis=0)
    shp = Shape(Geometry(coords, ["X"] * cn))
    rows = []
    for lab in refs:
        try:
            rows.append((float(shp.measure(lab, central_atom=0)), lab))
        except Exception:
            pass
    rows.sort()
    return rows

def uniq_metals(at, m_idx):
    ds = spglib.get_symmetry_dataset(
        (np.asarray(at.get_cell()), at.get_scaled_positions(), np.asarray(at.numbers)),
        symprec=1e-2)
    eq = ds.equivalent_atoms
    seen, reps = {}, []
    for mi in m_idx:
        k = int(eq[mi])
        if k not in seen:
            seen[k] = mi
            reps.append(mi)
    return reps, ds.international, len(reps)

def fmt_block(rows, n_show=None):
    r = rows if n_show is None else rows[:n_show]
    return "\n".join(f"      {lab:9s} {desc(lab):24s} S = {s:7.4f}{'  <--' if lab==rows[0][1] else ''}"
                     for s, lab in r)

def report(name, at, m_idx, jobs):
    """jobs: list of (include_I, cn, label). Each unique metal -> each job."""
    print(f"\n{'='*78}\n{name}   ({len(at)} atoms)")
    reps, spg, nu = uniq_metals(at, m_idx)
    print(f"  metals: {len(m_idx)} total, {nu} crystallographically unique (spglib: {spg})")
    for mi in reps:
        pos = at.get_positions()
        sym = at[mi].symbol
        for include_I, cn, label in jobs:
            dn = donors_of(at, mi, include_I)
            elem = [d[2] for d in dn]
            # enforce the intended CN: if donor count != cn, still report but flag
            tag = "" if len(elem) == cn else f"  [found {len(elem)} donors]"
            print(f"\n  >> {sym} (atom {mi})  [{label}]  donors={elem}  CN={len(elem)}{tag}")
            if len(elem) < 3 or len(elem) > 9:
                print("      CN out of SHAPE range [3..9]; skip"); continue
            mpos = at.get_positions()[mi]
            coords = np.array([mpos + dn[k][3] for k in range(len(dn))])
            rows = cshm_rows(coords, len(elem))
            print(fmt_block(rows))

def main():
    # ---- pristine (3 data blocks in one CIF: Y, Gd, Tb) ----
    pris = read(f"{DIR}/1907532-1999761.cif", index=":")
    metal_order = []
    for at in pris:
        ms = [i for i, s in enumerate(at.get_chemical_symbols()) if s in METALS]
        metal_order.append(at[ms[0]].symbol if ms else "?")
    # map block -> metal by counting: identify each block's metal
    for at in pris:
        at.pbc = True
        ms = [i for i, s in enumerate(at.get_chemical_symbols()) if s in METALS]
        el = at[ms[0]].symbol
        report(f"PRISTINE [{el}]", at, ms, [(False, 8, "6 framework O + 2 solvent O")])

    # ---- stripped (3 CIFs) ----
    for el, fn in [("Tb", "06_2019[Tb][pts]3[ASR]1.cif"),
                   ("Gd", "08_2020[Gd][pts]3[ASR]1.cif"),
                   ("Y",  "13_2020[Y][pts]3[ASR]1.cif")]:
        at = read(f"{DIR}/{fn}"); at.pbc = True
        ms = [i for i, s in enumerate(at.get_chemical_symbols()) if s in METALS]
        report(f"STRIPPED [{el}]", at, ms, [(False, 6, "6 framework O (solvent removed)")])

    # ---- POSCAR (3; true P1, all metals, CN=6 AND CN=7) ----
    for el, fn in [("Gd", "POSCAR_final_Gd"), ("Tb", "POSCAR_final_Tb"), ("Y", "POSCAR_final_Y")]:
        at = read(f"{DIR}/{fn}"); at.pbc = True
        ms = [i for i, s in enumerate(at.get_chemical_symbols()) if s in METALS]
        report(f"POSCAR [{el}] (DFT, I-loaded)", at, ms,
               [(False, 6, "framework only: 6 O"),
                (True,  7, "6 O + inner iodine")])

if __name__ == "__main__":
    main()
