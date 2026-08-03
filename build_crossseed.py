#!/usr/bin/env python3
"""Cross-seed Y and Tb with Gd's optimized uniform_chain iodine geometry.

Motivation: Gd's optimized chain has 2 very short I-O contacts (2.77 A) that Y/Tb lack;
Gd's electronic term is ~24 kJ/mol more negative. Test whether Y/Tb were trapped in a
shallower basin by seeding them with Gd's I arrangement.

Method: affine-register Gd's framework onto each target (Y, Tb) framework via least-squares
(M @ r_Gd + t ~ r_target, atom correspondence by nearest same-species), then apply M,t to
Gd's optimized I atoms -> target cell.  Framework of the target is frozen (its own POSCAR_final).
"""
import os, numpy as np, ase.io
import build_inputs as bi
BASE = "/workspace/dimitry/calc/qe/i2_mof"
SRC = "08_Gd/uniform_chain/POSCAR_final"   # Gd optimized chain (framework + I)
TARGETS = [("13_Y", "Y"), ("06_Tb", "Tb")]


def framework_and_I(path):
    a = ase.io.read(path)
    sy = np.array(a.get_chemical_symbols())
    fw = a[sy != "I"].copy()
    Iatoms = a[sy == "I"].copy()
    return fw, Iatoms


def fingerprint(atoms, n=12):
    """For each atom, sorted min-image distances to its n nearest neighbours (all atoms).
    Origin/rotation/cell-setting invariant -> robust correspondence."""
    cell = np.array(atoms.cell); inv = np.linalg.inv(cell); pos = atoms.positions
    m = len(pos)
    fp = np.zeros((m, n))
    for i in range(m):
        d = pos - pos[i]
        d = d - (cell.T @ np.round(inv @ d.T)).T
        ds = np.sort(np.linalg.norm(d, axis=1))[1:n+1]   # exclude self
        fp[i] = ds
    return fp


def correspond(src_fw, tgt_fw, src_metal, tgt_metal):
    """Match each src framework atom to the tgt atom of matching species (metals->'M')
    with the closest distance fingerprint. Returns tgt index per src atom."""
    ss = np.array([x if x != src_metal else "M" for x in src_fw.get_chemical_symbols()])
    ts = np.array([x if x != tgt_metal else "M" for x in tgt_fw.get_chemical_symbols()])
    fps = fingerprint(src_fw)
    fpt = fingerprint(tgt_fw)
    idx = np.zeros(len(ss), dtype=int)
    for k in range(len(ss)):
        cands = np.where(ts == ss[k])[0]
        err = np.sum((fpt[cands] - fps[k]) ** 2, axis=1)
        idx[k] = cands[int(np.argmin(err))]
    if len(np.unique(idx)) != len(idx):
        # resolve duplicates: keep best, remap losers to next-best unused
        used = set(); order = np.argsort([np.min(np.sum((fpt[np.where(ts==ss[k])[0]]-fps[k])**2,axis=1)) for k in range(len(ss))])
        for k in order:
            cands = np.where(ts == ss[k])[0]
            err = np.sum((fpt[cands] - fps[k]) ** 2, axis=1)
            for j in cands[np.argsort(err)]:
                if j not in used:
                    idx[k] = j; used.add(j); break
    return idx


def affine_map(src_fw, tgt_fw, idx):
    """Least-squares M(3x3),t so that tgt_pos ~ M@src_pos + t."""
    X = src_fw.positions
    Yp = tgt_fw.positions[idx]
    n = len(X)
    A = np.zeros((3 * n, 12))
    b = np.zeros(3 * n)
    for k in range(n):
        x = X[k]
        A[3*k:3*k+3, 0:9] = np.kron(np.eye(3), x)
        A[3*k:3*k+3, 9:12] = np.eye(3)
        b[3*k:3*k+3] = Yp[k]
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    M = sol[:9].reshape(3, 3, order="C").T   # row/col convention
    t = sol[9:12]
    pred = X @ M.T + t
    rms = np.sqrt(np.mean(np.sum((pred - Yp) ** 2, axis=1)))
    return M, t, rms


def main():
    os.chdir(BASE)
    src_fw, src_I = framework_and_I(SRC)
    print(f"Gd source: {len(src_fw)} fw atoms, {len(src_I)} I")
    for tgt, metal in TARGETS:
        tgt_fw, _ = framework_and_I(f"{tgt}/POSCAR_final")
        idx = correspond(src_fw, tgt_fw, "Gd", metal)
        M, t, rms = affine_map(src_fw, tgt_fw, idx)
        print(f"\n{tgt}: affine residual RMS = {rms*1000:.1f} mA (over {len(src_fw)} fw atoms)")
        # map Gd I -> target Cartesian, wrap into target cell
        cell = np.array(tgt_fw.cell)
        Icart = src_I.positions @ M.T + t
        inv = np.linalg.inv(cell)
        Ifrac = Icart @ inv.T
        Ifrac = np.mod(Ifrac, 1.0)
        # build target framework + mapped I
        new = tgt_fw.copy()
        new += ase.Atoms("I8", cell=cell, scaled_positions=Ifrac)
        outdir = f"{tgt}/from_Gd"
        os.makedirs(outdir, exist_ok=True)
        seed = f"{outdir}/POSCAR_seed"
        ase.io.write(seed, new, format="vasp", direct=True, vasp5=True)
        # clash check
        pos = new.positions; sy = np.array(new.get_chemical_symbols()); inv = np.linalg.inv(np.array(new.cell))
        def mic(i, j): d = pos[i]-pos[j]; d = d-np.array(new.cell).T@np.round(inv@d); return np.linalg.norm(d)
        Ii = [i for i, x in enumerate(sy) if x == "I"]; fw = [i for i, x in enumerate(sy) if x != "I"]
        minII = min(mic(i, j) for i in Ii for j in Ii if i > j)
        minIf = min(mic(i, j) for i in Ii for j in fw)
        print(f"  {tgt}/from_Gd: min I-I={minII:.2f} A, min I-framework={minIf:.2f} A")
        bi.build(seed, outdir, f"{tgt}_from_Gd")


if __name__ == "__main__":
    main()
