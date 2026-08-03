#!/usr/bin/env python3
"""Extract the final relaxed geometry from a Quantum ESPRESSO relax .out and write
a VASP POSCAR (direct/crystal coords) into the same directory.

QE prints "Begin final coordinates" only when BFGS converges. When the relax hits the
nstep cap (as most I2@MOF vdW relaxes do), that block is absent, so we fall back to the
LAST "ATOMIC_POSITIONS" block in the file -- the geometry of the final ionic step.

Cell is read from the matching .in (CELL_PARAMETERS); the framework is frozen so the cell
equals the original POSCAR cell.

Usage:  python qe_out_to_poscar.py <system>/loaded.out
"""
import re, sys
from pathlib import Path


def last_atomic_positions(text):
    lines = text.splitlines()
    # find the last "ATOMIC_POSITIONS (crystal)" header
    idx = max((i for i, ln in enumerate(lines) if ln.strip().startswith("ATOMIC_POSITIONS")), default=-1)
    if idx < 0:
        raise ValueError("No ATOMIC_POSITIONS block found")
    pos = []
    for line in lines[idx + 1:]:
        parts = line.split()
        if len(parts) < 4:
            if pos:
                break  # blank/short line after atoms -> end of block
            continue
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            break  # non-atom line (e.g. "Writing config-only ...") -> end of block
        sym = parts[0]
        if not sym.isalpha() or len(sym) > 2:
            break
        pos.append((sym, x, y, z))
    if not pos:
        raise ValueError("No atoms parsed after last ATOMIC_POSITIONS header")
    return pos


def lattice_from_input(in_file):
    text = in_file.read_text()
    m = re.search(r"CELL_PARAMETERS\s+\w+\s*\n(.*?)\nATOMIC_POSITIONS", text, re.DOTALL)
    rows = [list(map(float, ln.split())) for ln in m.group(1).strip().splitlines()]
    return rows


def write_poscar(pos, lattice, title, path):
    order, counts = [], {}
    for s, *_ in pos:
        if s not in counts:
            order.append(s); counts[s] = 0
        counts[s] += 1
    lines = [title, "1.0"]
    for r in lattice:
        lines.append(f"  {r[0]: .10f}  {r[1]: .10f}  {r[2]: .10f}")
    lines.append("  ".join(f"{s:>4}" for s in order))
    lines.append("  ".join(f"{counts[s]:>4}" for s in order))
    lines.append("Direct")
    for s, x, y, z in pos:
        lines.append(f"  {x:18.12f}  {y:18.12f}  {z:18.12f}")
    path.write_text("\n".join(lines) + "\n")


def main():
    out = Path(sys.argv[1])
    inp = out.with_suffix(".in")
    pos = last_atomic_positions(out.read_text())
    lat = lattice_from_input(inp)
    dest = out.parent / "POSCAR_final"
    write_poscar(pos, lat, f"final relaxed geometry from {out.name}", dest)
    print(f"{out.parent.name}: wrote {dest} ({len(pos)} atoms, {len(set(s for s, *_ in pos))} species)")


if __name__ == "__main__":
    main()
