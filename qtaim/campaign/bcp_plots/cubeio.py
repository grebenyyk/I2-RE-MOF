"""Cube + structure IO for the I2@MOF QTAIM contour plots.

Parses QE/critic2 cube headers (cell vectors, atoms), converts to
fractional / Angstrom, and provides MIC distance + neighbor helpers.
Coordinates in cube files are in bohr; cell vectors are voxel*count.
"""
import numpy as np

BOHR = 0.5291772109  # bohr -> Angstrom
EL = {"C": 6, "H": 1, "O": 8, "Y": 39, "Tb": 65, "Gd": 64, "I": 53, "N": 7, "Cd": 48}
ELZ = {v: k for k, v in EL.items()}


def read_cube_header(path):
    """Return cell_A (3x3 rows = a1,a2,a3 in Angstrom), atoms list of (sym, cart_A)."""
    with open(path) as f:
        lines = f.readlines()
    # first two lines are comments
    nat = int(lines[2].split()[0])
    cell = np.zeros((3, 3))
    nvox = []
    for i in range(3):
        p = lines[3 + i].split()
        nvox.append(int(p[0]))
        v = np.array([float(p[1]), float(p[2]), float(p[3])])  # bohr
        cell[i] = v * nvox[i] * BOHR
    atoms = []
    for k in range(nat):
        p = lines[6 + k].split()
        z = int(round(float(p[0])))
        r = np.array([float(p[2]), float(p[3]), float(p[4])]) * BOHR  # bohr->A
        atoms.append((ELZ.get(z, "?"), r))
    return cell, atoms, nvox


def read_cube_grid(path):
    """Read full cube: return cell_A, atoms, nvox, grid(n1,n2,n3)."""
    import os
    size = os.path.getsize(path)
    with open(path) as f:
        lines = f.readlines()
    nat = int(lines[2].split()[0])
    cell = np.zeros((3, 3)); nvox = []
    for i in range(3):
        p = lines[3 + i].split(); nvox.append(int(p[0]))
        v = np.array([float(p[1]), float(p[2]), float(p[3])])
        cell[i] = v * nvox[i] * BOHR
    atoms = []
    for k in range(nat):
        p = lines[6 + k].split()
        z = int(round(float(p[0])))
        r = np.array([float(p[2]), float(p[3]), float(p[4])]) * BOHR
        atoms.append((ELZ.get(z, "?"), r))
    # data starts at line 6+nat; join remainder and parse fast
    data_txt = "".join(lines[6 + nat:])
    data = np.fromstring(data_txt, sep=" ")
    grid = data.reshape(nvox[2], nvox[1], nvox[0]).transpose(2, 1, 0)  # ->(n1,n2,n3)
    return cell, atoms, nvox, grid


def interp_frac(grid, f):
    """Trilinear interpolation of a periodic cube grid at fractional coord f (len-3)."""
    n = np.array(grid.shape)
    g = f * n
    i0 = np.floor(g).astype(int) % n
    i1 = (i0 + 1) % n
    t = g - np.floor(g)
    c000 = grid[i0[0], i0[1], i0[2]]
    c100 = grid[i1[0], i0[1], i0[2]]
    c010 = grid[i0[0], i1[1], i0[2]]
    c110 = grid[i1[0], i1[1], i0[2]]
    c001 = grid[i0[0], i0[1], i1[2]]
    c101 = grid[i1[0], i0[1], i1[2]]
    c011 = grid[i0[0], i1[1], i1[2]]
    c111 = grid[i1[0], i1[1], i1[2]]
    c00 = c000 * (1 - t[0]) + c100 * t[0]
    c10 = c010 * (1 - t[0]) + c110 * t[0]
    c01 = c001 * (1 - t[0]) + c101 * t[0]
    c11 = c011 * (1 - t[0]) + c111 * t[0]
    c0 = c00 * (1 - t[1]) + c10 * t[1]
    c1 = c01 * (1 - t[1]) + c11 * t[1]
    return c0 * (1 - t[2]) + c1 * t[2]


def frac_from_cart(cell, cart):
    """cart (Angstrom) -> fractional using cell rows a1,a2,a3 (so cart = f . cell_rows)."""
    return np.linalg.solve(cell.T, cart)


def cart_from_frac(cell, f):
    return cell.T @ f


def mic_vec_frac(fc1, fc2, cell):
    """Minimum-image fractional displacement, returned as Angstrom Cartesian vector."""
    d = fc1 - fc2
    d = d - np.round(d)
    return cell.T @ d


def mic_dist(c1, c2, cell):
    """MIC distance between two Cartesian (Angstrom) points."""
    f1 = frac_from_cart(cell, c1)
    f2 = frac_from_cart(cell, c2)
    v = mic_vec_frac(f1, f2, cell)
    return np.linalg.norm(v)


def find_chain(struct_dir, cubename="val.cube", bond_max=3.6, neigh=6):
    """Print I geometry: I atoms (frac, cart-A), I-I pairs < bond_max,
    and each I's nearest O / metal neighbors. Returns dict."""
    import os
    path = os.path.join(struct_dir, cubename)
    cell, atoms, _nvox = read_cube_header(path)
    I_idx = [i for i, (s, _) in enumerate(atoms) if s == "I"]
    frac = {i: frac_from_cart(cell, atoms[i][1]) for i in range(len(atoms))}
    print(f"\n=== {struct_dir}  ({len(atoms)} atoms, {len(I_idx)} I) ===")
    print(f"cell (A): a={np.linalg.norm(cell[0]):.3f} b={np.linalg.norm(cell[1]):.3f} "
          f"c={np.linalg.norm(cell[2]):.3f}")
    print(f"\nI atoms:")
    for k, i in enumerate(I_idx, 1):
        f = frac[i]
        print(f"  I{k:2d} (at#{i+1:3d})  frac=({f[0]: .4f},{f[1]: .4f},{f[2]: .4f})  "
              f"cart=({atoms[i][1][0]:7.3f},{atoms[i][1][1]:7.3f},{atoms[i][1][2]:7.3f}) A")
    print(f"\nI-I pairs < {bond_max} A (MIC):")
    pairs = []
    for a in range(len(I_idx)):
        for b in range(a + 1, len(I_idx)):
            ia, ib = I_idx[a], I_idx[b]
            d = mic_dist(atoms[ia][1], atoms[ib][1], cell)
            if d < bond_max:
                pairs.append((d, a + 1, b + 1, ia, ib))
    for d, ka, kb, ia, ib in sorted(pairs):
        print(f"  I{ka}-I{kb}  d={d:.3f} A")
    print(f"\nnearest framework neighbors of each I (within {neigh} A):")
    metals = {"Y", "Tb", "Gd", "Cd"}
    for k, i in enumerate(I_idx, 1):
        nbrs = []
        for j, (s, r) in enumerate(atoms):
            if j == i:
                continue
            d = mic_dist(atoms[i][1], r, cell)
            if d < neigh and (s in metals or s == "O"):
                nbrs.append((d, s, j + 1))
        nbrs.sort()
        brief = ", ".join(f"{s}#{j}={d:.2f}" for d, s, j in nbrs[:4])
        print(f"  I{k}: {brief}")
    return {"cell": cell, "atoms": atoms, "frac": frac, "I_idx": I_idx, "pairs": pairs}


if __name__ == "__main__":
    import sys
    base = "/workspace/dimitry/calc/qe/i2_mof/qtaim/campaign"
    for s in (sys.argv[1:] or ["Y_str", "Y_act", "Y_actc"]):
        find_chain(f"{base}/{s}")
