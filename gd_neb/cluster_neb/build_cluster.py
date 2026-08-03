#!/usr/bin/env python3
"""Build a finite cluster around the ch1 iodines for NWChem NEB.

ch1 (molecular atom order, 1-based): I1, I3, I4, I5.
In pw_1.in the I atoms appear in file order; indices 1..8 correspond to
atoms I1..I8 of the molecular numbering. So ch1 = (I1,I3,I4,I5) = the
1st, 3rd, 4th, 5th iodine in pw_1.in.

Steps:
 1. Read pw_1.in (molecular) and pw_7.in (chain) with ASE -> Atoms (cart).
    Framework atoms (C/Gd/H/O) are identical in both; only iodines differ.
 2. Use the framework from pw_1 (molecular) — frozen, identical for both
    endpoints.
 3. Cart coords for ch1 iodines: molecular from pw_1; chain from pw_7.
 4. Pick framework atoms within CUT Ang of ANY ch1 iodine (use the union
    of molecular+chain iodine positions so the wall encloses both states).
 5. Cap dangling C-C and C-O bonds at the cut with H (1.09 A C-H, 0.96 A O-H)
    along the broken-bond direction.
 6. Write xyz files (mol + chain) and a constraints/frozen-atom index map.
"""
import numpy as np
from ase.io import read, write
from ase.neighborlist import NeighborList, natural_cutoffs

CUT = 3.8           # Ang, framework atom kept if within CUT of any ch1 I
                    # (tight wall: just the nearest carboxylate O + aromatic C
                    #  that directly confine the iodines; cheaper SCF)
CH1 = [1, 3, 4, 5]  # 1-based I indices within the I-block (1..8)
DROP_GD = True      # remove Gd from QM region (avoids f7 open-shell + basis gap)
GD_CHARGE = 3.0     # Gd(III)
USE_POINTCHARGES = False  # with neutral protonated carboxylates, drop Gd PCs
                          # (they'd create an artificial field near neutral COOH)

PW1 = "/workspace/dimitry/calc/qe/i2_mof/gd_neb/pw_1.in"
PW7 = "/workspace/dimitry/calc/qe/i2_mof/gd_neb/pw_7.in"
OUT = "/workspace/dimitry/calc/qe/i2_mof/gd_neb/cluster_neb"


def parse_qe_in(path):
    """Parse ATOMIC_POSITIONS {crystal|angstrom} from a pw.x input."""
    cell = []
    species = []
    coords = []
    pos_block = None
    cell_block = False
    reading_species = False
    n_per = {}
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        s = line.strip()
        low = s.lower()
        if low.startswith("cell_parameters"):
            # next lines until blank/keyword
            cell_block = True
            i += 1
            while i < len(lines):
                t = lines[i].split()
                if len(t) >= 3:
                    try:
                        cell.append([float(x) for x in t[:3]])
                        i += 1
                        continue
                    except ValueError:
                        break
                break
            continue
        if low.startswith("atomic_positions"):
            pos_block = "crystal" if "crystal" in low else ("angstrom" if "angstrom" in low else "crystal")
            i += 1
            while i < len(lines):
                t = lines[i].split()
                if len(t) >= 4 and t[0].isalpha():
                    try:
                        species.append(t[0])
                        coords.append([float(t[1]), float(t[2]), float(t[3])])
                        i += 1
                        continue
                    except ValueError:
                        break
                break
            continue
        i += 1
    cell = np.array(cell)
    coords = np.array(coords)
    return cell, species, coords, pos_block


def to_atoms(path):
    """Read pw_*.in via ASE (it understands pw.scf in format='espresso-in')."""
    return read(path, format="espresso-in")


mol = to_atoms(PW1)
chn = to_atoms(PW7)
print("mol atoms:", len(mol), "symbols:", sorted(set(mol.get_chemical_symbols())))
print("chn atoms:", len(chn), "symbols:", sorted(set(chn.get_chemical_symbols())))

# Identify iodine indices (in pw_*.in the species are grouped C,Gd,H,O,I).
I_idx = [i for i, s in enumerate(mol.get_chemical_symbols()) if s == "I"]
print("I indices (0-based) in mol:", I_idx)
ch1_idx = [I_idx[k - 1] for k in CH1]
print("ch1 indices (0-based):", ch1_idx)

# Framework atoms = everything that is not I.
fw_idx = [i for i, s in enumerate(mol.get_chemical_symbols()) if s != "I"]
mol_i_ch1 = mol.positions[ch1_idx]
chn_i_ch1 = chn.positions[ch1_idx]

# Build the "envelope" iodine positions: union of mol & chain ch1 coords,
# so the wall is big enough to enclose both endpoints.
i_env = np.vstack([mol_i_ch1, chn_i_ch1])

# PBC-safe distance: cell is small (~9x9x17 with tilt). Use min-image.
cell = mol.cell.array
inv = np.linalg.inv(cell.T)


def mic_diff(p, q):
    d = p - q
    f = d @ inv
    f -= np.round(f)
    return f @ cell


# Select framework atoms within CUT of ANY env iodine (mic).
keep_fw = []
for j in fw_idx:
    p = mol.positions[j]
    dmin = min(np.linalg.norm(mic_diff(p, ip)) for ip in i_env)
    if dmin <= CUT:
        keep_fw.append(j)
print(f"framework atoms within {CUT} A of any ch1 I: {len(keep_fw)}")

# Optionally drop Gd from QM region (captured as point charges instead).
gd_positions = []
if DROP_GD:
    keep_fw = [j for j in keep_fw if mol.get_chemical_symbols()[j] != "Gd"]
    # collect ALL Gd in the unit cell for point-charge embedding (4 of them)
    gd_positions = [mol.positions[j] for j, s in enumerate(mol.get_chemical_symbols()) if s == "Gd"]
    print(f"DROPPED Gd from QM; {len(gd_positions)} Gd point charges added at +3")

# Compose the cluster: framework atoms (frozen) + 4 ch1 iodines (mobile).
# Order: framework first, then iodines.
from ase import Atoms

fw_symbols = [mol.get_chemical_symbols()[j] for j in keep_fw]
fw_pos = mol.positions[keep_fw]

cluster_fw = Atoms(symbols=fw_symbols, positions=fw_pos)


def cap_dangling(atoms):
    """Cap only dangling C bonds at the cluster boundary with H.
    Carboxylate O that coordinated a (now-point-charge) Gd is left as-is
    (the +3 point charge provides electrostatic compensation), so we do
    NOT cap O — that would wrongly convert carboxylate O- into hydroxyl.
    """
    sym = atoms.get_chemical_symbols()
    pos = atoms.get_positions()
    cutoffs = natural_cutoffs(atoms, mult=1.2)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)
    caps = []
    cap_len = {"C": 1.09, "O": 0.96}
    for i, s in enumerate(sym):
        if s not in ("C", "O"):
            continue
        neigh, _ = nl.get_neighbors(i)
        n = len(neigh)
        if s == "C":
            expected = 3
        else:  # O
            expected = 2  # carboxylate/ether O has 2 heavy-atom neighbors;
                          # dangling carboxylate O (lost Gd) has 1 -> cap with H (->COOH)
        missing = expected - n
        if missing <= 0:
            continue
        if n == 0:
            dirs = [np.array([1.0, 0.0, 0.0])]
        else:
            vsum = np.zeros(3)
            for k in neigh:
                vsum += pos[k] - pos[i]
            vsum /= np.linalg.norm(vsum) + 1e-9
            base = -vsum
            dirs = []
            for m in range(missing):
                if missing == 1:
                    d = base
                else:
                    perp = np.array([0.0, 0.0, 1.0]) if abs(base[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
                    perp = np.cross(base, perp)
                    perp /= np.linalg.norm(perp) + 1e-9
                    ang = (m - (missing - 1) / 2.0) * 0.6
                    d = np.cos(ang) * base + np.sin(ang) * perp
                dirs.append(d)
        for d in dirs:
            d = d / (np.linalg.norm(d) + 1e-9)
            hp = pos[i] + cap_len[s] * d
            caps.append((hp, i))
    new_atoms = atoms.copy()
    for hp, pi in caps:
        new_atoms.append("H")
        new_atoms.positions[-1] = hp
    return new_atoms, len(caps)


# Cap on framework only.
cluster_fw_capped, n_caps = cap_dangling(cluster_fw)
print(f"added {n_caps} capping H atoms")

# Now build the two endpoint clusters: capped frozen framework + 4 ch1 iodines.
mol_cluster = cluster_fw_capped.copy()
mol_cluster += Atoms(symbols=["I"] * 4, positions=mol_i_ch1)

chn_cluster = cluster_fw_capped.copy()
chn_cluster += Atoms(symbols=["I"] * 4, positions=chn_i_ch1)

print("mol_cluster: natoms =", len(mol_cluster), "symbols:", sorted(set(mol_cluster.get_chemical_symbols())))
print("chn_cluster: natoms =", len(chn_cluster))

# Frozen set: all framework+capping atoms are frozen; only last 4 (iodines) move.
n_total = len(mol_cluster)
frozen = list(range(n_total - 4))  # 0-based frozen indices
mobile = list(range(n_total - 4, n_total))

# Sanity: I-I distances in each endpoint (mic not needed, cluster is finite)
def pair_d(atoms):
    p = atoms.positions[-4:]
    out = {}
    for a in range(4):
        for b in range(a + 1, 4):
            d = np.linalg.norm(p[a] - p[b])
            if d < 4.0:
                out[f"I{a+1}-I{b+1}"] = round(d, 3)
    return out

print("mol short I-I:", pair_d(mol_cluster))
print("chn short I-I:", pair_d(chn_cluster))

write(f"{OUT}/cluster_mol.xyz", mol_cluster, format="xyz")
write(f"{OUT}/cluster_chn.xyz", chn_cluster, format="xyz")

# Write point-charge file for NWChem (BQ charges at Gd positions).
if gd_positions and USE_POINTCHARGES:
    with open(f"{OUT}/pointcharges.q", "w") as qf:
        qf.write(f"# Gd(III) point charges (electrostatic embedding for dropped Gd)\n")
        for p in gd_positions:
            qf.write(f"{GD_CHARGE:.1f} {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}\n")
    print(f"wrote pointcharges.q with {len(gd_positions)} Gd point charges")

# Save frozen/mobile indices for downstream NWChem input builder.
np.savetxt(f"{OUT}/frozen_indices.txt", frozen, fmt="%d")
with open(f"{OUT}/cluster_info.txt", "w") as f:
    f.write(f"natoms {n_total}\n")
    f.write(f"frozen {len(frozen)} mobile {len(mobile)}\n")
    f.write(f"mobile_indices_0based {' '.join(map(str, mobile))}\n")
    f.write(f"ch1_i_positions_mol:\n")
    for p in mol_i_ch1:
        f.write(f"  {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
    f.write(f"ch1_i_positions_chn:\n")
    for p in chn_i_ch1:
        f.write(f"  {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")

print("Wrote cluster_mol.xyz, cluster_chn.xyz, cluster_info.txt, frozen_indices.txt")
