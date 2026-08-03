#!/usr/bin/env python3
"""ASE NEB driver for the molecular->chain I4 transformation.

Uses NWChem (PBE+D3, def2-SVP on framework, def2-SVP+ECP on I) as the
force engine. Framework atoms are held fixed by ASE (constraints), so
NWChem sees an unconstrained cluster but ASE only propagates the 4
mobile iodines along the band.

Strategy:
  - Read optimized endpoints (opt_mol_*.xyz, opt_chn_*.xyz) if present,
    else fall back to the unrelaxed cluster endpoints.
  - 7 images, IDPP initial interpolation (more robust than linear for
    this kind of bond rearrangement).
  - NEB with spring constant 0.1 eV/A^2, climb-image after a few iters.
  - Each image -> its own subdir so NWChem scratch doesn't collide.
"""
import os, sys, shutil, subprocess, time
import numpy as np
from ase.io import read, write
from ase.neb import NEB, IDPP
from ase.optimize import BFGS, LBFGS
from ase import Atoms

OUT = "/workspace/dimitry/calc/qe/i2_mof/gd_neb/cluster_neb"
NIMAGES = 7       # interior images (excluding endpoints) -> total = NIMAGES+2
SPRING = 0.1      # eV/A^2
FMAX = 0.10       # eV/A  NEB convergence
CI = True         # climbing image


def find_nwchem():
    # use the nwchem conda env's nwchem
    return "/workspace/miniconda3/envs/nwchem/bin/nwchem"


def env_setup():
    e = os.environ.copy()
    e["NWCHEM_BASIS_LIBRARY"] = "/workspace/miniconda3/envs/nwchem/share/nwchem/libraries/"
    e["NWCHEM_NWPW_LIBRARY"] = "/workspace/miniconda3/envs/nwchem/share/nwchem/libraryps/"
    e["OMP_NUM_THREADS"] = "1"
    e["OMPI_ALLOW_RUN_AS_ROOT"] = "1"
    e["OMPI_ALLOW_RUN_AS_ROOT_CONFIRM"] = "1"
    e["PATH"] = "/workspace/miniconda3/envs/nwchem/bin:" + e["PATH"]
    return e


def nwchem_template(xc="xpbe96 cpbe96", disp="disp vdw 3", i_basis="def2-svp"):
    return f"""
memory stack 600 mb heap 1200 mb global 2400 mb
charge 0

dft
  mult 1
  xc {xc}
  iterations 100
  convergence energy 1.0e-07
  convergence density 1.0e-06
  convergence gradient 1.0e-04
  grid fine
  mulliken
  {disp}
end

basis
  C   library def2-svp
  H   library def2-svp
  O   library def2-svp
  I   library {i_basis}
end
ecp
  I   library def2-ecp
end

geometry noautosym nocenter units angstroms
<<GEOM>>
end

task dft energy
"""


def write_nwchem_in(atoms, path, template):
    lines = []
    for s, p in zip(atoms.get_chemical_symbols(), atoms.positions):
        lines.append(f"  {s:2s} {p[0]:14.8f} {p[1]:14.8f} {p[2]:14.8f}")
    body = template.replace("<<GEOM>>", "\n".join(lines))
    with open(path, "w") as f:
        f.write(body)


def parse_nwchem_energy(out_path):
    """Get the last 'Total DFT energy' value from an NWChem output."""
    e = None
    if not os.path.exists(out_path):
        return None
    with open(out_path) as f:
        for line in f:
            if "Total DFT energy" in line:
                try:
                    e = float(line.split()[-1])
                except Exception:
                    pass
    return e


def parse_nwchem_forces(out_path, natoms):
    """Parse the final NWChem DFT energy gradient block into (F, E).
    Looks for 'DFT ENERGY GRADIENTS' block with columns:
       atom ... x y z (the gradient). NWChem prints gradients; ASE wants forces = -gradient.
    """
    if not os.path.exists(out_path):
        return None, None
    F = np.zeros((natoms, 3))
    E = parse_nwchem_energy(out_path)
    in_grad = False
    rows = []
    with open(out_path) as f:
        for line in f:
            if "DFT ENERGY GRADIENTS" in line:
                in_grad = True
                rows = []
                continue
            if in_grad:
                if line.strip().startswith("-----") or "------------" in line:
                    if rows:
                        in_grad = False
                        break
                    continue
                tok = line.split()
                # format: idx sym nuc charge x y z (gradient components)
                if len(tok) >= 7 and tok[1].isalpha():
                    try:
                        gx, gy, gz = float(tok[-3]), float(tok[-2]), float(tok[-1])
                        rows.append((int(tok[0]), gx, gy, gz))
                    except ValueError:
                        pass
                elif line.strip() == "" and rows:
                    in_grad = False
                    break
    for idx, gx, gy, gz in rows:
        if 1 <= idx <= natoms:
            F[idx - 1] = [-gx, -gy, -gz]
    if not rows:
        return None, E
    return F, E


def run_nwchem_sp(atoms, workdir, label, template, env, nwchem_exe, natoms, nproc=16, reuse_db=False):
    """Run an NWChem single point in workdir; return (energy, forces)."""
    os.makedirs(workdir, exist_ok=True)
    in_path = os.path.join(workdir, f"{label}.nw")
    out_path = os.path.join(workdir, f"{label}.out")
    write_nwchem_in(atoms, in_path, template)
    # NWChem leaves movecs/db files; clean to avoid reuse surprises unless requested.
    if not reuse_db:
        for ext in (".db", ".movecs", ".gridpts.00"):
            p = os.path.join(workdir, label + ext)
            if os.path.exists(p):
                os.remove(p)
    cmd = ["mpirun", "--allow-run-as-root", "--oversubscribe", "-np", str(nproc), nwchem_exe, in_path]
    t0 = time.time()
    try:
        with open(out_path, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT, env=env, cwd=workdir, timeout=1800)
    except subprocess.TimeoutExpired:
        return None, None
    dt = time.time() - t0
    E = parse_nwchem_energy(out_path)
    F, _ = parse_nwchem_forces(out_path, natoms)
    return E, F, dt


class NWChemCalc:
    def __init__(self, template, env, nwchem_exe, root, nproc=16):
        self.template = template
        self.env = env
        self.nwchem_exe = nwchem_exe
        self.root = root
        self.nproc = nproc
        self.counter = 0

    def get_energy_and_forces(self, atoms, label):
        natoms = len(atoms)
        wd = os.path.join(self.root, f"img_{label}")
        E, F, dt = run_nwchem_sp(
            atoms, wd, label, self.template, self.env, self.nwchem_exe, natoms, self.nproc
        )
        return E, F, dt


def manual_neb_step(images_pos, endpoints_E, calc, mobile_mask, spring, ci_image_idx=None, k_info=None):
    """One NEB step is too custom; we'll use ASE NEB instead below.
    Kept for reference; not used.
    """
    raise NotImplementedError


def main():
    # Endpoints: prefer optimized xyz if present
    mol_path = os.path.join(OUT, "opt_mol_pbe_svp_opt.xyz")
    chn_path = os.path.join(OUT, "opt_chn_pbe_svp_opt.xyz")
    if not os.path.exists(mol_path):
        mol_path = os.path.join(OUT, "cluster_mol.xyz")
    if not os.path.exists(chn_path):
        chn_path = os.path.join(OUT, "cluster_chn.xyz")
    print(f"endpoints: mol={mol_path}\n          chn={chn_path}")
    mol = read(mol_path)
    chn = read(chn_path)
    natoms = len(mol)
    mobile = list(range(natoms - 4, natoms))
    print(f"natoms={natoms} mobile_iodines(0-based)={mobile}")

    # Build initial band with IDPP interpolation on the mobile atoms only.
    # ASE NEB wants a list of atoms; endpoints fixed by setting their constraints.
    images = [mol]
    for i in range(1, NIMAGES + 1):
        img = mol.copy()
        images.append(img)
    images.append(chn)

    # IDPP init: only interpolate the mobile atoms; copy framework from mol.
    # Easiest: use ASE's IDPP on the full images but fix the framework via constraints.
    # Actually we just need linear interp of the 4 iodines; framework identical.
    mol_p = mol.positions
    chn_p = chn.positions
    for k in range(1, NIMAGES + 1):
        frac = k / (NIMAGES + 1)
        pos = mol_p.copy()
        for m in mobile:
            pos[m] = mol_p[m] * (1 - frac) + chn_p[m] * frac
        images[k].set_positions(pos)

    # Pre-compute endpoint energies
    template = nwchem_template()
    env = env_setup()
    nwchem_exe = find_nwchem()
    calc = NWChemCalc(template, env, nwchem_exe, os.path.join(OUT, "neb_scratch"), nproc=16)

    print("=== computing endpoint energies ===")
    E0, F0, _ = calc.get_energy_and_forces(mol, "endpoint_mol")
    print(f"E_mol = {E0}")
    EN, FN, _ = calc.get_energy_and_forces(chn, "endpoint_chn")
    print(f"E_chn = {EN}")

    if E0 is None or EN is None:
        print("ERROR: endpoint SP failed; aborting.")
        return

    # Manual NEB loop (simple elastic band + steepest descent on mobile atoms).
    # This is a lightweight NEB; convergence is slow but robust and we control it.
    spring_k = SPRING  # eV/A^2
    lr = 0.05           # step size (A) per iteration on real forces (tuned for safety)
    lr_spring = 1.0     # spring relaxation step (much faster; springs are harmonic)
    max_iter = 40

    # Current positions of mobile atoms in each interior image
    interior_pos = [images[k].positions[mobile].copy() for k in range(1, NIMAGES + 1)]
    # Endpoints (fixed)
    left_mobile = mol_p[mobile].copy()
    right_mobile = chn_p[mobile].copy()

    energies = [E0]
    history = []

    for it in range(max_iter):
        print(f"\n--- NEB iteration {it+1}/{max_iter} ---")
        # 1) compute forces on each interior image (parallel would be nice; serial for now)
        image_forces = []
        image_energies = []
        for k in range(1, NIMAGES + 1):
            img = mol.copy()
            full_pos = mol_p.copy()
            full_pos[mobile] = interior_pos[k - 1]
            img.set_positions(full_pos)
            E, F, dt = calc.get_energy_and_forces(img, f"neb_img{k}")
            image_energies.append(E)
            image_forces.append(F)
            print(f"  img{k}: E={E}  ({dt:.0f}s)" if E is not None else f"  img{k}: FAILED")
        if any(E is None for E in image_energies):
            print("ERROR: an image failed; stopping.")
            break

        # 2) find climbing image (highest energy)
        ci = int(np.argmax(image_energies))
        print(f"  climbing image = img{ci+1} (E={image_energies[ci]:.6f})")

        # 3) build tangent + spring force on each interior image's mobile atoms
        all_mobile = [left_mobile] + interior_pos + [right_mobile]
        new_interior = [p.copy() for p in interior_pos]
        max_f_perp = 0.0
        for k in range(1, NIMAGES + 1):
            pos = interior_pos[k - 1]
            prev = all_mobile[k - 1]
            nxt = all_mobile[k + 1]
            # tangent (normalized), using energy-weighted choice (Henkelman)
            d_prev = pos - prev
            d_next = nxt - pos
            n_prev = np.linalg.norm(d_prev)
            n_next = np.linalg.norm(d_next)
            E_prev = energies[-1] if k == 1 else image_energies[k - 2]
            E_cur = image_energies[k - 1]
            E_next = EN if k == NIMAGES else image_energies[k]
            if E_next > E_cur > E_prev:
                tau = d_next / (n_next + 1e-12)
            elif E_prev > E_cur > E_next:
                tau = -d_prev / (n_prev + 1e-12)
            else:
                tau = (d_next / (n_next + 1e-12) - d_prev / (n_prev + 1e-12))
                tau /= np.linalg.norm(tau) + 1e-12
            # real force on mobile atoms
            Freal = image_forces[k - 1][mobile]  # (4,3)
            # split into perpendicular and parallel
            Fpar = np.sum(Freal * tau, axis=1, keepdims=True) * tau
            Fperp = Freal - Fpar
            # spring force (parallel only)
            Fspring = spring_k * (d_next - d_prev)[:, None] * np.eye(3)[None]  # not quite; do componentwise:
            Fspring = spring_k * ((nxt - pos) - (pos - prev)) / (np.linalg.norm(nxt - pos) + np.linalg.norm(pos - prev) + 1e-12)[:, None]
            # zero spring perpendicular component
            spar = np.sum(Fspring * tau, axis=1, keepdims=True) * tau
            Fspring = spar
            # total force
            if CI and (k - 1) == ci:
                # climbing image: invert parallel real force, no spring
                Ftotal = Fperp - Fpar
            else:
                Ftotal = Fperp + Fspring
            max_f_perp = max(max_f_perp, np.max(np.linalg.norm(Ftotal, axis=1)))
            # update (steepest descent)
            new_interior[k - 1] = pos + lr * Ftotal
        # commit
        for k in range(1, NIMAGES + 1):
            interior_pos[k - 1] = new_interior[k - 1]
        energies_iter = [E0] + image_energies + [EN]
        energies = energies_iter  # for next-iter tangent energy lookup
        history.append(energies_iter)
        barrier = max(image_energies) - E0
        # relative to lower endpoint
        Eref = min(E0, EN)
        rel = [e - Eref for e in energies_iter]
        print(f"  barrier vs mol = {barrier*27.2114:.4f} eV = {barrier*2625.5:.1f} kJ/mol")
        print(f"  max|F_neb| on mobile = {max_f_perp:.3f} eV/A")
        print(f"  E profile (kJ/mol vs min endpoint): " + " ".join(f"{r*2625.5:6.1f}" for r in rel))
        # save band
        with open(os.path.join(OUT, f"neb_band_iter{it+1:02d}.xyz"), "w") as bf:
            for k in range(NIMAGES + 2):
                img = mol.copy()
                full_pos = mol_p.copy()
                full_pos[mobile] = (left_mobile if k == 0 else (right_mobile if k == NIMAGES + 1 else interior_pos[k - 1]))
                img.set_positions(full_pos)
                bf.write(f"{natoms}\nimage {k} iter {it+1} E={energies_iter[k]}\n")
                for s, p in zip(img.get_chemical_symbols(), img.positions):
                    bf.write(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
        with open(os.path.join(OUT, "neb_profile.txt"), "a") as pf:
            pf.write(f"iter {it+1} maxF {max_f_perp:.4f} " + " ".join(f"{r*2625.5:.2f}" for r in rel) + "\n")
        if max_f_perp < FMAX:
            print(f"  CONVERGED at Fmax={max_f_perp:.3f} < {FMAX}")
            break

    # final report
    print("\n=== FINAL NEB profile ===")
    final = history[-1]
    Eref = min(E0, EN)
    rel = [(e - Eref) * 2625.5 for e in final]
    print("kJ/mol vs min endpoint:", " ".join(f"{r:7.2f}" for r in rel))
    print("eV     vs min endpoint:", " ".join(f"{r/2625.5*27.2114:7.3f}" for r in rel))
    barrier_mol = (max(final) - E0) * 2625.5
    barrier_chn = (max(final) - EN) * 2625.5
    print(f"barrier (TS - mol) = {barrier_mol:.1f} kJ/mol = {barrier_mol/96.485:.3f} eV")
    print(f"barrier (TS - chn) = {barrier_chn:.1f} kJ/mol = {barrier_chn/96.485:.3f} eV")


if __name__ == "__main__":
    main()
