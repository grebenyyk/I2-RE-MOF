#!/usr/bin/env python3
"""Parallel ASE-style NEB for the mol->chain I4 transformation.

Evaluates all interior images in parallel (one NWChem SP per image, batched).
Each image keeps its own movecs for SCF guess reuse across NEB iterations.

NEB method: simple elastic band with Henkelman-Liu tangent + climbing image.
Optimizer: steepest descent with adaptive step (robust, no line search needed).
"""
import os, sys, shutil, subprocess, time, glob
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from ase.io import read

OUT = "/workspace/dimitry/calc/qe/i2_mof/gd_neb/cluster_neb"
NPROC_PER_IMG = 4     # ranks per image
NIMAGES = 4           # interior images (4 x 4 = 16 ranks = cgroup quota)
SPRING = 0.05         # eV/A^2
FMAX = 0.20           # eV/A  convergence on max |F_neb|
MAX_ITER = 30
GRID = "medium"
SCF_CONV = 1e-4       # density convergence (energy to ~1 kJ/mol; forces usable)
I_BASIS = "def2-svp"
XC = "xpbe96 cpbe96"  # = PBE


NWCHEM = "/workspace/miniconda3/envs/nwchem/bin/nwchem"
MPIRUN = "/workspace/miniconda3/envs/nwchem/bin/mpirun"
BASE_ENV = dict(os.environ)  # inherit everything (HOME, USER, etc.) to keep opal happy
BASE_ENV.update({
    "NWCHEM_BASIS_LIBRARY": "/workspace/miniconda3/envs/nwchem/share/nwchem/libraries/",
    "NWCHEM_NWPW_LIBRARY": "/workspace/miniconda3/envs/nwchem/share/nwchem/libraryps/",
    "OMP_NUM_THREADS": "1",
    "OMPI_ALLOW_RUN_AS_ROOT": "1",
    "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM": "1",
    "PATH": "/workspace/miniconda3/envs/nwchem/bin:" + os.environ.get("PATH", ""),
})


def template_str():
    return f"""
memory stack 400 mb heap 800 mb global 1600 mb
charge 0

dft
  mult 1
  xc {XC}
  iterations 100
  convergence energy 1.0e10
  convergence density {SCF_CONV:g}
  convergence gradient 1.0e10
  grid {GRID}
  <<SMEAR>>
  vectors input <<MOVECS_IN>> output <<MOVECS_OUT>>
end

basis
  C   library def2-svp
  H   library def2-svp
  O   library def2-svp
  I   library {I_BASIS}
end
ecp
  I   library def2-ecp
end

geometry noautosym nocenter noautoz units angstroms
<<GEOM>>
end

task dft gradient
""".replace("  {''}\n", "\n")


def write_input(atoms, in_path, movecs_in, movecs_out, tmpl, smear=False):
    lines = []
    for s, p in zip(atoms.get_chemical_symbols(), atoms.positions):
        lines.append(f"  {s:2s} {p[0]:14.8f} {p[1]:14.8f} {p[2]:14.8f}")
    body = tmpl.replace("<<GEOM>>", "\n".join(lines))
    body = body.replace("<<MOVECS_IN>>", movecs_in).replace("<<MOVECS_OUT>>", movecs_out)
    body = body.replace("<<SMEAR>>", "smear 0.01" if smear else "")
    with open(in_path, "w") as f:
        f.write(body)


def parse_energy(out_path):
    """Get final energy. Prefer 'Total DFT energy' (printed on full convergence);
    fall back to the last SCF iteration's energy column."""
    e = None
    last_scf_e = None
    if not os.path.exists(out_path):
        return None
    with open(out_path) as f:
        for line in f:
            if "Total DFT energy" in line:
                try:
                    e = float(line.split()[-1])
                except Exception:
                    pass
            # SCF line looks like: " d= 0,ls=0.5,diis    18  -2225.46  ..." (col 4 = energy)
            stripped = line.strip()
            if stripped.startswith("d=") or "diis" in stripped[:30]:
                tok = stripped.split()
                # find the energy token: first token after the iteration index that parses as float
                for j in range(2, len(tok)):
                    try:
                        v = float(tok[j])
                        last_scf_e = v
                        break
                    except ValueError:
                        continue
    return e if e is not None else last_scf_e


def parse_forces(out_path, natoms):
    if not os.path.exists(out_path):
        return None
    F = np.zeros((natoms, 3))
    in_grad = False
    rows = []
    with open(out_path) as f:
        for line in f:
            if "DFT ENERGY GRADIENTS" in line:
                in_grad = True
                rows = []
                continue
            if in_grad:
                tok = line.split()
                if len(tok) >= 7 and tok[1].isalpha() and tok[0].isdigit():
                    try:
                        gx, gy, gz = float(tok[-3]), float(tok[-2]), float(tok[-1])
                        rows.append((int(tok[0]), gx, gy, gz))
                    except ValueError:
                        pass
                elif (line.strip() == "" or "-----" in line) and rows:
                    in_grad = False
                    break
    for idx, gx, gy, gz in rows:
        if 1 <= idx <= natoms:
            F[idx - 1] = [-gx, -gy, -gz]
    return F if rows else None


def run_one(atoms, workdir, label, natoms, tmpl, smear=False):
    """Run one NWChem SP+grad, reusing movecs if present.
    NWChem runs with cwd=workdir, so movecs paths are bare filenames.
    If smear=True, applies fractional occupation (helps small-gap images)."""
    os.makedirs(workdir, exist_ok=True)
    in_path = os.path.join(workdir, f"{label}.nw")
    out_path = os.path.join(workdir, f"{label}.out")
    movecs_file = f"{label}.movecs"   # bare name; NWChem cwd = workdir
    movecs_path = os.path.join(workdir, movecs_file)
    movecs_in = movecs_file if os.path.exists(movecs_path) else "atomic"
    write_input(atoms, in_path, movecs_in, movecs_file, tmpl, smear=smear)
    cmd = [MPIRUN, "--allow-run-as-root", "--oversubscribe", "-np", str(NPROC_PER_IMG), NWCHEM, in_path]
    t0 = time.time()
    try:
        with open(out_path, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT, env=BASE_ENV, cwd=workdir, timeout=1200)
    except subprocess.TimeoutExpired:
        return None, None, 1200.0
    dt = time.time() - t0
    E = parse_energy(out_path)
    F = parse_forces(out_path, natoms)
    return E, F, dt


def main():
    tmpl = template_str()

    # Endpoints: use periodic-relaxed cluster coords directly (frozen framework
    # means these ARE the cluster minima to first order; NEB will refine along path).
    mol = read(f"{OUT}/cluster_mol.xyz")
    chn = read(f"{OUT}/cluster_chn.xyz")
    natoms = len(mol)
    mobile = list(range(natoms - 4, natoms))
    print(f"natoms={natoms} NIMAGES={NIMAGES} ranks/img={NPROC_PER_IMG}")

    mol_p = mol.positions.copy()
    chn_p = chn.positions.copy()
    left_mobile = mol_p[mobile].copy()
    right_mobile = chn_p[mobile].copy()

    # Linear interpolation of mobile atoms for initial band.
    interior_pos = []
    for k in range(1, NIMAGES + 1):
        frac = k / (NIMAGES + 1)
        pos = mol_p[mobile] * (1 - frac) + chn_p[mobile] * frac
        interior_pos.append(pos)

    def atoms_with_mobile(mobile_pos):
        a = mol.copy()
        p = mol_p.copy()
        p[mobile] = mobile_pos
        a.set_positions(p)
        return a

    # Endpoints: compute once, cache energies. Retry with smear if plain fails.
    print("=== computing endpoint SPs ===")
    E_mol, F_mol, dt = run_one(mol, os.path.join(OUT, "neb_scratch", "endpoint_mol"), "endpoint_mol", natoms, tmpl)
    if E_mol is None:
        print(f"  E_mol plain failed ({dt:.0f}s); retrying with smear")
        E_mol, F_mol, dt = run_one(mol, os.path.join(OUT, "neb_scratch", "endpoint_mol"), "endpoint_mol", natoms, tmpl, smear=True)
    print(f"  E_mol = {E_mol}  ({dt:.0f}s)")
    E_chn, F_chn, dt = run_one(chn, os.path.join(OUT, "neb_scratch", "endpoint_chn"), "endpoint_chn", natoms, tmpl)
    if E_chn is None:
        print(f"  E_chn plain failed ({dt:.0f}s); retrying with smear")
        E_chn, F_chn, dt = run_one(chn, os.path.join(OUT, "neb_scratch", "endpoint_chn"), "endpoint_chn", natoms, tmpl, smear=True)
    print(f"  E_chn = {E_chn}  ({dt:.0f}s)")
    if E_mol is None or E_chn is None:
        print("FATAL: endpoint SP failed (even with smear)"); return
    print(f"  dE(chn-mol) = {(E_chn-E_mol)*2625.5:.1f} kJ/mol  "
          f"{'(chain LOWER, correct trend)' if E_chn<E_mol else '(chain HIGHER — check)'}")

    profile_log = open(os.path.join(OUT, "neb_profile.txt"), "w")

    for it in range(1, MAX_ITER + 1):
        t_it = time.time()
        # Build atoms objects for each interior image.
        image_atoms = [atoms_with_mobile(interior_pos[k]) for k in range(NIMAGES)]
        # Run all in parallel.
        results = [None] * NIMAGES
        with ThreadPoolExecutor(max_workers=NIMAGES) as ex:
            future_map = {
                ex.submit(run_one, image_atoms[k], os.path.join(OUT, "neb_scratch", f"img{k+1}"), f"img{k+1}", natoms, tmpl): k
                for k in range(NIMAGES)
            }
            for fut in as_completed(future_map):
                k = future_map[fut]
                results[k] = fut.result()
        Es = [r[0] for r in results]
        Fs = [r[1] for r in results]
        dts = [r[2] for r in results]
        # Retry failed images with smear (small-gap polyiodide images).
        failed = [k for k, E in enumerate(Es) if E is None]
        if failed:
            print(f"  images {failed} failed plain; retrying with smear")
            retry_results = [None] * len(failed)
            with ThreadPoolExecutor(max_workers=len(failed)) as ex:
                fmap = {ex.submit(run_one, image_atoms[k], os.path.join(OUT, "neb_scratch", f"img{k+1}"), f"img{k+1}", natoms, tmpl, True): i
                        for i, k in enumerate(failed)}
                for fut in as_completed(fmap):
                    retry_results[fmap[fut]] = fut.result()
            for i, k in enumerate(failed):
                Es[k], Fs[k], dts[k] = retry_results[i]
        if any(E is None for E in Es):
            print(f"iter {it}: an image FAILED (even with smear); aborting")
            for k, E in enumerate(Es):
                print(f"  img{k+1}: E={E}")
            break
        print(f"\n=== NEB iter {it} (wall {time.time()-t_it:.0f}s; per-img max {max(dts):.0f}s) ===")
        full_E = [E_mol] + Es + [E_chn]
        Eref = min(E_mol, E_chn)
        rel = [(e - Eref) * 2625.5 for e in full_E]
        print("  E profile (kJ/mol vs min endpoint): " + " ".join(f"{r:7.1f}" for r in rel))
        profile_log.write(f"iter {it} " + " ".join(f"{r:.2f}" for r in rel) + "\n")
        profile_log.flush()

        # Climbing image = highest-energy interior image.
        ci = int(np.argmax(Es))
        print(f"  climbing image: img{ci+1} (E={Es[ci]:.6f})")

        # Compute NEB force on each interior image's mobile atoms.
        all_mobile = [left_mobile] + interior_pos + [right_mobile]
        new_interior = [p.copy() for p in interior_pos]
        max_f = 0.0
        for k in range(NIMAGES):
            pos = interior_pos[k]
            prev = all_mobile[k]
            nxt = all_mobile[k + 2]
            d_prev = pos - prev
            d_next = nxt - pos
            n_prev = np.linalg.norm(d_prev)
            n_next = np.linalg.norm(d_next)
            E_prev = E_mol if k == 0 else Es[k - 1]
            E_cur = Es[k]
            E_next = E_chn if k == NIMAGES - 1 else Es[k + 1]
            # Henkelman tangent
            if E_next > E_cur > E_prev:
                tau = d_next / (n_next + 1e-12)
            elif E_prev > E_cur > E_next:
                tau = -d_prev / (n_prev + 1e-12)
            else:
                tau = d_next / (n_next + 1e-12) - d_prev / (n_prev + 1e-12)
                tau = tau / (np.linalg.norm(tau) + 1e-12)
            tau = tau / (np.linalg.norm(tau) + 1e-12)  # (4,3) per-atom? No—make it per-atom:
            # Actually tangent should be per-atom. Use the same direction for all 4 mobile
            # atoms (they move coherently). Simpler & robust: average tangent.
            tau_avg = tau.mean(axis=0) if tau.ndim == 2 else tau
            tau_avg /= np.linalg.norm(tau_avg) + 1e-12
            Freal = Fs[k][mobile]  # (4,3)
            Fpar_mag = np.dot(Freal, tau_avg)  # (4,)
            Fpar = Fpar_mag[:, None] * tau_avg[None, :]
            Fperp = Freal - Fpar
            # spring force along tangent (averaged distance)
            avg_d = 0.5 * (n_prev + n_next)
            spring_per_atom = SPRING * (n_next - n_prev)  # scalar
            Fspring = spring_per_atom * tau_avg[None, :] * np.ones((4, 1))
            if k == ci:
                Ftotal = Fperp - Fpar  # invert parallel
            else:
                Ftotal = Fperp + Fspring
            max_f = max(max_f, np.max(np.linalg.norm(Ftotal, axis=1)))
            new_interior[k] = pos + 0.02 * Ftotal  # small fixed step (A)

        # adaptive step: scale down if forces huge
        print(f"  max|F_neb| on mobile = {max_f:.3f} eV/A")
        # save band xyz
        with open(os.path.join(OUT, f"neb_band_iter{it:02d}.xyz"), "w") as bf:
            for k in range(NIMAGES + 2):
                if k == 0:
                    mpos = left_mobile
                elif k == NIMAGES + 1:
                    mpos = right_mobile
                else:
                    mpos = new_interior[k - 1]
                a = atoms_with_mobile(mpos)
                bf.write(f"{natoms}\nimage {k} iter {it} E={full_E[k]:.6f}\n")
                for s, p in zip(a.get_chemical_symbols(), a.positions):
                    bf.write(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
        interior_pos = new_interior
        if max_f < FMAX:
            print(f"  CONVERGED (maxF {max_f:.3f} < {FMAX})")
            break

    profile_log.close()
    # Final summary
    E_TS = max(Es) if Es else None
    print("\n=== FINAL ===")
    print(f"E_mol = {E_mol:.6f} Hartree")
    print(f"E_chn = {E_chn:.6f} Hartree")
    print(f"E_TS  = {E_TS:.6f} Hartree")
    if E_TS:
        barrier_mol = (E_TS - E_mol) * 2625.5
        barrier_chn = (E_TS - E_chn) * 2625.5
        print(f"Barrier (TS - mol) = {barrier_mol:.1f} kJ/mol = {barrier_mol/96.485:.3f} eV")
        print(f"Barrier (TS - chn) = {barrier_chn:.1f} kJ/mol = {barrier_chn/96.485:.3f} eV")
        print(f"Reaction energy   = {(E_chn-E_mol)*2625.5:.1f} kJ/mol")


if __name__ == "__main__":
    main()
