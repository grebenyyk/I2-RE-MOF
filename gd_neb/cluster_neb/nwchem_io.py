#!/usr/bin/env python3
"""Build NWChem input files for the cluster NEB.

Single point + geometry optimization (with frozen framework) + helper to
write a series of NEB image inputs.

Method:
  - DFT PBE0-D3(BJ) [fallback PBE+D3 if needed]
  - def2-SVP on C/H/O, def2-ECP+def2-SVP on Gd (4f in core), def2-ECP+def2-TZVP on I
  - RI not used (NWChem RI for DFT needs aux basis setup; skip for robustness)
  - SCF: tight (1e-6), 100 iterations, restricted closed-shell
  - Geometry opt: only the 4 iodines move (use constraints block).
"""
import sys, os
from ase.io import read

OUT = "/workspace/dimitry/calc/qe/i2_mof/gd_neb/cluster_neb"


def geom_block(atoms, name="cluster", cartesian_opt=True):
    """Return NWChem geometry block with Cartesian coords.
    noautoz forces Cartesian coords (required when using 'fix atom' constraints)."""
    kw = "noautosym nocenter noautoz" if cartesian_opt else "noautosym nocenter"
    lines = [f"geometry {kw} units angstroms"]
    lines.append(f"  # {name}")
    sym = atoms.get_chemical_symbols()
    pos = atoms.positions
    for s, p in zip(sym, pos):
        lines.append(f"  {s:2s} {p[0]:14.8f} {p[1]:14.8f} {p[2]:14.8f}")
    lines.append("end")
    return "\n".join(lines)


def basis_block(i_basis="def2-tzvp"):
    """def2-SVP on framework (C/H/O), {i_basis}+ECP on iodine."""
    return f"""
basis
  C   library def2-svp
  H   library def2-svp
  O   library def2-svp
  I   library {i_basis}
end
ecp
  I   library def2-ecp
end
"""


def dft_header(scf_conv=1e-6, maxiter=100, xc="pbe0", disp="disp vdw 3", grid="medium"):
    return f"""
memory stack 600 mb heap 1200 mb global 2400 mb
charge 0

dft
  mult 1
  xc {xc}
  iterations {maxiter}
  convergence energy 1.0e-06
  convergence density {scf_conv:g}
  convergence gradient 1.0e-04
  grid {grid}
  mulliken
  {disp}
end
"""


def write_sp(atoms, name, xc="pbe0", disp="disp vdw 3", i_basis="def2-tzvp"):
    body = (
        dft_header(xc=xc, disp=disp)
        + basis_block(i_basis=i_basis)
        + "\n"
        + geom_block(atoms, name)
        + "\n\ntask dft energy\n"
    )
    path = f"{OUT}/{name}.nw"
    with open(path, "w") as f:
        f.write(body)
    return path


def write_opt(atoms, name, mobile_mask, xc="pbe0", disp="disp vdw 3", maxiter=80, i_basis="def2-tzvp"):
    """Geometry opt with framework frozen. mobile_mask = list of 0-based indices that may move."""
    n = len(atoms)
    mobile = set(mobile_mask)
    cons = ["constraints"]
    for i in range(n):
        if i not in mobile:
            cons.append(f"  fix atom {i+1}")  # NWChem atoms are 1-based
    cons.append("end")
    body = (
        dft_header(xc=xc, disp=disp)
        + basis_block(i_basis=i_basis)
        + "\n"
        + geom_block(atoms, name)
        + "\n\n"
        + "\n".join(cons)
        + f"""

driver
  maxiter {maxiter}
  xyz {name}_opt
end

task dft optimize
"""
    )
    path = f"{OUT}/{name}.nw"
    with open(path, "w") as f:
        f.write(body)
    return path


if __name__ == "__main__":
    mol = read(f"{OUT}/cluster_mol.xyz")
    chn = read(f"{OUT}/cluster_chn.xyz")
    mobile = list(range(len(mol) - 4, len(mol)))

    # Production NEB method: PBE+D3 (matches periodic functional), def2-SVP on I for speed.
    p = write_sp(mol, "sp_mol_pbe_svp", xc="xpbe96 cpbe96", disp="disp vdw 3", i_basis="def2-svp")
    print("wrote", p)
    p = write_sp(chn, "sp_chn_pbe_svp", xc="xpbe96 cpbe96", disp="disp vdw 3", i_basis="def2-svp")
    print("wrote", p)

    # Endpoint optimizations (frozen framework, relax 4 iodines), SVP-on-I
    p = write_opt(mol, "opt_mol_pbe_svp", mobile, xc="xpbe96 cpbe96", disp="disp vdw 3", i_basis="def2-svp", maxiter=100)
    print("wrote", p)
    p = write_opt(chn, "opt_chn_pbe_svp", mobile, xc="xpbe96 cpbe96", disp="disp vdw 3", i_basis="def2-svp", maxiter=100)
    print("wrote", p)

    # Also PBE0+D3 single points on the endpoints (with TZVP on I) for a hybrid-func cross-check
    p = write_sp(mol, "sp_mol_pbe0", xc="pbe0", disp="disp vdw 3", i_basis="def2-tzvp")
    print("wrote", p)
    p = write_sp(chn, "sp_chn_pbe0", xc="pbe0", disp="disp vdw 3", i_basis="def2-tzvp")
    print("wrote", p)

    # Keep the original PBE+D3 TZVP-on-I smoke test too
    p = write_sp(mol, "sp_mol_pbe", xc="xpbe96 cpbe96", disp="disp vdw 3", i_basis="def2-tzvp")
    print("wrote", p)
