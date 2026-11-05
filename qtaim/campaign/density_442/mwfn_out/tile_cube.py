"""Tile a QE cube along one lattice axis to extend the periodic density.

Standard mode: nrep copies along `axis`, coverage [0, nrep*A], atoms at +k*A.
  --neg : instead cover [-A, +A] (2 copies) by shifting the origin to -A and
          appending the -A image cell. The ORIGINAL atoms keep their indices
          (listed first), so atom-referenced planes/labels are unchanged. This
          fills a blank on the negative side of the plot window.

Cube convention: z (n3) fastest -> reshape (n1,n2,n3); x slowest.
"""
import sys
import numpy as np


def tile(cube_in, cube_out, nrep=2, axis=0, neg=False):
    lines = open(cube_in).readlines()
    title1, title2 = lines[0], lines[1]
    p3 = lines[2].split()
    nat = int(p3[0]); origin = np.array([float(x) for x in p3[1:4]])

    def vec(ln):
        p = ln.split(); return int(p[0]), np.array([float(p[1]), float(p[2]), float(p[3])])
    n1, h1 = vec(lines[3]); n2, h2 = vec(lines[4]); n3, h3 = vec(lines[5])
    ns = [n1, n2, n3]; hs = [h1, h2, h3]
    A_axis = ns[axis] * hs[axis]
    atom_lines = lines[6:6 + nat]

    data = np.array("".join(lines[6 + nat:]).split(), dtype=float)
    assert data.size == n1 * n2 * n3, f"grid {data.size} != {n1*n2*n3}"
    arr = data.reshape(n1, n2, n3)

    if neg:                                  # cover [-A,+A]: 2 copies, origin -A
        shifts = [np.zeros(3), -A_axis]      # originals (idx preserved) + -A image
        big = np.concatenate([arr, arr], axis=axis)
        out_origin = origin - A_axis
        factor = 2
    else:
        shifts = [k * A_axis for k in range(nrep)]
        big = np.concatenate([arr] * nrep, axis=axis)
        out_origin = origin
        factor = nrep

    with open(cube_out, "w") as o:
        o.write(title1); o.write(title2)
        o.write(f"{nat*len(shifts):5d} {out_origin[0]:12.6f} {out_origin[1]:12.6f} {out_origin[2]:12.6f}\n")
        for ax in range(3):
            nn = ns[ax] * factor if ax == axis else ns[ax]
            hh = hs[ax]
            o.write(f"{nn:5d} {hh[0]:12.6f} {hh[1]:12.6f} {hh[2]:12.6f}\n")
        for sh in shifts:
            for ln in atom_lines:
                p = ln.split()
                x = np.array([float(p[2]), float(p[3]), float(p[4])]) + sh
                o.write(f"{int(p[0]):5d} {float(p[1]):12.6f} {x[0]:12.6f} {x[1]:12.6f} {x[2]:12.6f}\n")
        flat = big.reshape(-1)
        for i in range(0, flat.size, 6):
            o.write(" ".join(f"{v:.5E}" for v in flat[i:i + 6]) + "\n")
    g = [ns[a] * (factor if a == axis else 1) for a in range(3)]
    print(f"wrote {cube_out}: axis={axis} neg={neg} -> {nat*len(shifts)} atoms, grid {g}, "
          f"origin {out_origin.round(3)}")


if __name__ == "__main__":
    a = sys.argv
    ci, co = a[1], a[2]
    nrep, axis, neg = 2, 0, False
    i = 3
    while i < len(a):
        if a[i].isdigit(): nrep = int(a[i])
        elif a[i] == "--axis": axis = int(a[i + 1]); i += 1
        elif a[i] == "--neg": neg = True
        i += 1
    tile(ci, co, nrep, axis, neg)
