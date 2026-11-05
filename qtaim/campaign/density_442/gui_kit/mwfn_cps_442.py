#!/usr/bin/env python3
"""Convert 442 critic2 outputs -> Multiwfn CPs.txt for the plane maps.

Row format (critic2 1.4.309; col2 is the point group, NOT an atom label):
  idx pg (3,sig) type fx fy fz mult name f |grad| lap
CP sources: bcp_raw_{tag}.out (auto scan) UNION refine_{tag}.out (seed-point
refinement that recovered the I-I BCPs the scan missed), deduped by MIC
proximity. Selection mirrors the 442 tables' classifier (nearest-2 nuclei,
MIC, from the unit-cell cube atoms):
  - BOND CPs only (dots); no nuclear / ring / cage
  - keep if either of the 2 nearest nuclei is iodine  (=> union of the II, IO,
    IC and MI tables)
  - drop if the OTHER nearest nucleus is H  (the C-H...I contacts the tables
    rejected as artifacts by path verification)
  - near-core filters as the tables: d2/d1 < 2 (kills the one-sided PAW I-4d
    shell artifacts) and d1 < 2.6 A (no floating CPs)
  - rho >= RHO_MIN a.u. (default 0.003: keeps the weak I...I cross-links
    0.0036-0.0071 and the I...O BCPs)
  - replicate CPs over the tiled supercell; per axis pass the integer image
    shifts matching tile_cube.py (originals first, indices unchanged), e.g.
    st: --shift0 0,1,2 --shift1 0,-1 --shift2 0,-1
"""
import re, sys
import numpy as np

BOHR = 0.5291772109
ROWRE = re.compile(
    r"^\s*(\d+)\s+\S+\s+\(\s*3\s*,\s*(-\d+)\s*\)\s+bond\s+"
    r"([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+\d+\s+\S+\s+([0-9.Ee+-]+)")


def cube_cell(cube):
    with open(cube) as f:
        f.readline(); f.readline()
        nat = int(f.readline().split()[0])
        v, zs, pos = [], [], []
        for _ in range(3):
            p = f.readline().split()
            v.append([int(p[0]) * float(p[1]), int(p[0]) * float(p[2]), int(p[0]) * float(p[3])])
        A = np.array(v)                                   # rows a1,a2,a3 (bohr)
        for _ in range(nat):
            p = f.readline().split()
            zs.append(int(float(p[0])))
            pos.append([float(p[2]), float(p[3]), float(p[4])])
    return A, np.array(pos), zs


def main(bcpouts, unit_cube, out, shifts, rho_min=0.003, ratio=2.0, max_d=2.6):
    A, pos, zs = cube_cell(unit_cube)
    Ain = np.linalg.inv(A)
    nucl = pos @ Ain
    cand = []
    for bcpout in bcpouts:
        for l in open(bcpout, errors='ignore'):
            m = ROWRE.match(l)
            if not m:
                continue
            cand.append((np.array([float(m.group(i)) for i in (3, 4, 5)]),
                         float(m.group(6))))
    uniq = []                                             # dedup < 0.15 A (MIC)
    for f, rho in cand:
        if not any(np.linalg.norm((f - g - np.round(f - g)) @ A * BOHR) < 0.15
                   for g, _ in uniq):
            uniq.append((f, rho))
    kept = []
    for f, rho in uniq:
        d = f - nucl; d -= np.round(d)
        dist = np.linalg.norm(d @ A * BOHR, axis=1)
        o = np.argsort(dist)
        d1, d2, k1, k2 = dist[o[0]], dist[o[1]], o[0], o[1]
        if zs[k1] == 1 or zs[k2] == 1:
            continue                                      # H contact (artifact)
        if zs[k1] != 53 and zs[k2] != 53:
            continue                                      # no iodine involved
        if d2 / d1 >= ratio or d1 > max_d:
            continue                                      # near-core / floating
        if rho < rho_min:
            continue
        kept.append((f, rho))
    rows = []                                             # replicate the tiling
    for s1 in shifts[0]:
        for s2 in shifts[1]:
            for s3 in shifts[2]:
                sh = s1 * A[0] + s2 * A[1] + s3 * A[2]    # bohr, tiled-cube frame
                for f, _ in kept:
                    rows.append(f @ A + sh)
    with open(out, 'w') as o:
        o.write(f'{len(rows):6d}\n')
        for i, c in enumerate(rows, 1):
            o.write(f'{i:6d} {c[0]:12.6f} {c[1]:12.6f} {c[2]:12.6f} {2:4d}\n')
    print(f'{out}: {len(kept)} iodine BCPs/cell (from {len(cand)} merged rows, '
          f'{len(uniq)} unique) x {len(rows) // max(len(kept), 1)} copies = {len(rows)} dots')


if __name__ == '__main__':
    a = sys.argv[1:]
    unit_cube, out = a[0], a[1]
    bcpouts, shifts = [], [None, None, None]
    i = 2
    while i < len(a):
        if a[i].startswith('--shift'):
            ax = int(a[i][7]); shifts[ax] = [int(x) for x in a[i + 1].split(',')]; i += 2
        else:
            bcpouts.append(a[i]); i += 1
    for k in range(3):
        if shifts[k] is None:
            shifts[k] = [0]
    main(bcpouts, unit_cube, out, shifts)
