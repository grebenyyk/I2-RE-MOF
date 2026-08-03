#!/usr/bin/env python3
"""
Build two UNIFORM channel-filling variants from each pts POSCAR_final (13_Y, 06_Tb, 08_Gd).

The relaxed (seeded) POSCAR_final packs iodine NON-uniformly: an extended polyiodide CHAIN
in one ab-layer (bridging that layer's two metals) and discrete I2 dimers (HERRINGBONE) in
the other ab-layer. The framework is P2_1/c, so the two layers are inversion partners.

Variants (test whether a uniform fill is deeper than the mixed minimum):
  uniform_chain     = chain-atoms + inversion(chain-atoms)         # chain in BOTH layers
  uniform_herring   = herring-atoms + inversion(herring-atoms)     # I2 in BOTH layers

Iodine classification is by connectivity (I-I < 3.6 A, min-image): the largest connected
component is the chain; the remaining dimer(s) are the herringbone. Inversion is the P2_1/c
op with R = -I, taken per-structure from spglib on the bare framework.

Writes, per variant: <system>/uniform_chain/{POSCAR_seed, loaded.in} and
                     <system>/uniform_herring/{POSCAR_seed, loaded.in}.
The framework is frozen (if_pos 0 0 0); only I relaxes (if_pos 1 1 1) — same template as the
campaign (david, smearing mv degauss 0.01, mixing 0.2, conv 1e-5).
"""
import os, numpy as np, ase.io, spglib
import build_inputs as bi

STRUCTS = [("13_Y", "Y"), ("06_Tb", "Tb"), ("08_Gd", "Gd")]
BASE = "/workspace/dimitry/calc/qe/i2_mof"
IICUT = 3.6   # I-I contact threshold (chain/dimer connectivity), A


def connected_components(Ii, pos, cell):
    inv = np.linalg.inv(cell)
    n = len(Ii)
    parent = list(range(n))
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b): parent[find(a)] = find(b)
    for a in range(n):
        for b in range(a + 1, n):
            d = pos[Ii[a]] - pos[Ii[b]]
            d = d - cell.T @ np.round(inv @ d)
            if np.linalg.norm(d) < IICUT:
                union(a, b)
    comps = {}
    for i in range(n):
        comps.setdefault(find(i), []).append(i)
    return list(comps.values())


def inversion_op(fw_atoms, symprec=1e-2):
    ds = spglib.get_symmetry_dataset(
        (np.array(fw_atoms.cell), fw_atoms.get_scaled_positions(), fw_atoms.get_atomic_numbers()),
        symprec=symprec)
    target = -np.eye(3, dtype=int)
    for R, t in zip(ds["rotations"], ds["translations"]):
        if np.array_equal(R, target):
            return t
    raise RuntimeError("no inversion in framework symmetry")


def build_variant(system, metal):
    atoms = ase.io.read(f"{system}/POSCAR_final")
    cell = np.array(atoms.cell)
    sy = atoms.get_chemical_symbols()
    f = atoms.get_scaled_positions(wrap=False)
    fw_idx = [i for i, s in enumerate(sy) if s != "I"]
    Ii = [i for i, s in enumerate(sy) if s == "I"]
    assert len(Ii) == 8, f"{system}: expected 8 I, got {len(Ii)}"

    comps = connected_components(Ii, atoms.positions, cell)
    comps.sort(key=len, reverse=True)
    chain_local = comps[0]                 # largest component = polyiodide chain
    herring_local = [j for c in comps[1:] for j in c]   # dimer(s) = herringbone
    chain_I = [Ii[k] for k in chain_local]
    herring_I = [Ii[k] for k in herring_local]
    print(f"\n[{system}] chain({len(chain_I)} I) largest-comp; herring({len(herring_I)} I) dimers")

    t_inv = inversion_op(atoms[fw_idx])
    def invert(idxset):
        out = []
        for i in idxset:
            p = -f[i] + t_inv              # x' = R@x + t with R=-I
            out.append(np.mod(p, 1.0))
        return np.array(out)

    for name, keep in [("uniform_chain", chain_I), ("uniform_herring", herring_I)]:
        keep_f = f[keep]
        inv_f = invert(keep)
        allI = np.vstack([keep_f, inv_f])
        assert len(allI) == 8
        # new atoms: framework + 8 I
        new = atoms[fw_idx].copy()
        ii = ase.Atoms("I8", cell=cell, scaled_positions=allI)
        new += ii
        outdir = f"{system}/{name}"
        os.makedirs(outdir, exist_ok=True)
        seed = f"{outdir}/POSCAR_seed"
        ase.io.write(seed, new, format="vasp", direct=True, vasp5=True)
        # sanity: no I-I clash, I count, chain/dimer preserved
        ipos = new.get_scaled_positions(wrap=False)[np.array(new.get_chemical_symbols()) == "I"]
        compsv = connected_components(list(range(8)),
                                      (np.array(new.cell).T @ ipos.T).T, cell)
        compv = sorted(len(c) for c in compsv)
        # generate relax input via the shared builder (framework frozen, I movable)
        bi.build(seed, outdir, f"{system}_{name}")
        print(f"  {name}: wrote {seed} + loaded.in | I-I components: {compv}")


if __name__ == "__main__":
    os.chdir(BASE)
    for sysname, metal in STRUCTS:
        build_variant(sysname, metal)
