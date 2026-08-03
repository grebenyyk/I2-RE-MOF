"""Convert critic2 bcp.out critical points -> Multiwfn CPs.txt

Load via Multiwfn main fn 2 -> -4 -> 5 so a plane map (main fn 4) overlays CP
dots (Multiwfn's own grid CP-search finds nothing on a cube).

bcp.out reports CPs in crystallographic (fractional) coords; CPs.txt wants Bohr
Cartesian. We read lattice vectors from the cube header and convert f@a.

CPs.txt format (verified via Multiwfn "Add a CP" -> "Save CPs"):
    <N>
    <idx> <x_bohr> <y_bohr> <z_bohr> <type>     # 1=NCP 2=BCP 3=ring 4=cage

FILTERS (curate the noise critic2 leaves near heavy-atom cores):
  ratio test   drop any CP whose 2 nearest nuclei are one-sided (d2/d1 >= RATIO);
               kills near-core spurious BCP/ring/cage ("one BCP looks wrong",
               "orange dots on atoms"). Real CPs sit between >=2 nuclei (d2/d1~1).
  ring/cage    additionally dropped if d1 < RC_FLOOR Angstrom (too close to a
               nucleus to be a real ring/cage centre).
  rho floor    drop CP with electron density rho < RHO_MIN (weak/long contacts,
               e.g. a 3.6 Angstrom I...I the chemist does not want drawn).
Nuclear CPs are always skipped (atoms are drawn anyway).
"""
import sys, re, numpy as np

BOHR = 0.5291772109
SIG2TYPE = {-3: 1, -1: 2, 1: 3, 3: 4}
ROWRE = re.compile(
    r"\(\s*3\s*,\s*([+-]?\d+)\s*\)\s+(?:nucleus|bond|ring|cage)\s+"
    r"([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+\d+\s+\S+\s+([0-9.Ee+-]+)")


def cube_cell(cube):
    with open(cube) as f:
        f.readline(); f.readline()
        nat = int(f.readline().split()[0])
        v = []
        for _ in range(3):
            p = f.readline().split()
            n = int(p[0]); v.append([n*float(p[1]), n*float(p[2]), n*float(p[3])])
        A = np.array(v)                       # rows a1,a2,a3 in Bohr
        atoms = np.array([[float(x) for x in f.readline().split()[2:5]]
                          for _ in range(nat)])  # Bohr Cartesian
    return A, atoms


def parse_bcpout(bcpout):
    """Yield (typecode, frac[3], rho) for every CP in the bcp.out CP table."""
    rows = []
    in_tbl = False
    for ln in open(bcpout, errors="ignore"):
        s = ln.strip()
        if s.startswith("# ncp") and "cryst" in s:
            in_tbl = True; continue
        if s.startswith("# ncp"):            # next table (connectivity) -> stop
            in_tbl = False
        if not in_tbl:
            continue
        m = ROWRE.search(s)
        if not m:
            continue
        sig = int(m.group(1))
        if sig not in SIG2TYPE:
            continue
        f = np.array([float(m.group(2)), float(m.group(3)), float(m.group(4))])
        rows.append((SIG2TYPE[sig], f, float(m.group(5))))
    return rows


def main(bcpout, cube, out="CPs.txt", keep_ring_cage=False,
         ratio=2.0, rho_min=0.0, rc_floor=1.0, tile_a1=False, tile_a2=0, a2_neg=False,
         max_d=2.6, tile_a3=0, a3_neg=False):
    A, atoms = cube_cell(cube)                            # lattice + atoms (bohr)
    Ain = np.linalg.inv(A)
    rows = parse_bcpout(bcpout)
    # use ALL cube atoms (full 80) as nuclei, not critic2's symmetry-reduced nuclear-CP list
    nucl = atoms @ Ain                                    # fractional nuclei

    def nearest2_ang(frac):
        d = frac - nucl
        d -= np.round(d)                                  # MIC wrap
        ds = np.sort(np.linalg.norm(d @ A * BOHR, axis=1))
        return ds[0], ds[1]

    kept = []
    for tc, f, rho in rows:
        if tc == 1:
            continue                                      # skip nuclear
        d1, d2 = nearest2_ang(f)
        if d2 / d1 >= ratio:                              # one-sided core artifact
            continue
        if d1 > max_d:                                    # floating (no atom within max_d)
            continue
        if tc in (3, 4):
            if not keep_ring_cage:
                continue
            if d1 < rc_floor:                             # too close to a nucleus
                continue
        if rho < rho_min:                                 # weak / long contact
            continue
        c = f @ A                                         # Bohr Cartesian (unit cell)
        kept.append((tc, c, rho))

    # duplicate CPs to match a tiled supercell along a1, a2, a3 (a2/a3 may be negative)
    reps = [int(tile_a1) if tile_a1 else 1,
            int(tile_a2) if tile_a2 else 1,
            int(tile_a3) if tile_a3 else 1]
    if any(r > 1 for r in reps):
        a1, a2, a3 = A[0], A[1], A[2]
        s2 = -1 if a2_neg else 1
        s3 = -1 if a3_neg else 1
        kept = [(tc, c + k1*a1 + k2*s2*a2 + k3*s3*a3, rho)
                for k1 in range(reps[0]) for k2 in range(reps[1]) for k3 in range(reps[2])
                for tc, c, rho in kept]

    with open(out, "w") as o:
        o.write(f"{len(kept):6d}\n")
        for i, (tc, c, _) in enumerate(kept, 1):
            o.write(f"{i:6d} {c[0]:12.6f} {c[1]:12.6f} {c[2]:12.6f} {tc:4d}\n")
    n = {}
    for tc, _, _ in kept:
        n[tc] = n.get(tc, 0) + 1
    print(f"wrote {out}: {len(kept)} CPs "
          f"(BCP={n.get(2,0)}, ring={n.get(3,0)}, cage={n.get(4,0)}); "
          f"ratio<{ratio}, rho>={rho_min}, rc_floor={rc_floor}A")


def _arg(name, default, argv):
    for a in argv:
        if a.startswith(name + "="):
            return type(default)(a.split("=", 1)[1])
    return default


def _tile(name, a):
    for x in a:
        if x == name:
            return 2
        if x.startswith(name + "="):
            return int(x.split("=", 1)[1])
    return 0


if __name__ == "__main__":
    a = sys.argv
    bcpout = a[1] if len(a) > 1 else "bcp.out"
    cube   = a[2] if len(a) > 2 else "val.cube"
    out    = a[3] if len(a) > 3 else "CPs.txt"
    main(bcpout, cube, out,
         keep_ring_cage="--with-ring-cage" in a,
         ratio=_arg("--ratio", 2.0, a),
         rho_min=_arg("--rho-min", 0.0, a),
         rc_floor=_arg("--rc-floor", 1.0, a),
         tile_a1=_tile("--tile-a1", a),
         tile_a2=_tile("--tile-a2", a),
         tile_a3=_tile("--tile-a3", a),
         a2_neg="--a2-neg" in a,
         a3_neg="--a3-neg" in a)
