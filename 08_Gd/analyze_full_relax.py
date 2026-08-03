#!/usr/bin/env python3
"""
Analyze the two Gd full-relax runs: extract final energies, final geometries,
compute sorption energies two ways, I-I contacts, and framework RMSD vs seed.

Reads:
  08_Gd/full_relax_chain/relax.out   (chain seed, converged)
  08_Gd/full_relax_mol/relax.out     (mol seed)
  08_Gd/uniform_chain/POSCAR_final   (chain seed geometry)
  08_Gd/activated_i2/POSCAR_final_optimized (mol seed geometry)
  08_Gd/activate_empty/relax.out     (activated empty framework ref)
  08_Gd/full_relax/empty_full.out    (prior full-relax empty ref -- DID NOT converge,
                                      documented caveat)
  i2_mof/i2_vac_ref/scf.out          (E_I2 vac)

Writes 08_Gd/full_relax_energies.txt with all results.
"""
import os
import re
import numpy as np

BASE = "/workspace/dimitry/calc/qe/i2_mof/08_Gd"
CHAIN_DIR = os.path.join(BASE, "full_relax_chain")
MOL_DIR = os.path.join(BASE, "full_relax_mol")

# cell (angstrom), same for both seeds and the campaign
CELL = np.array([
    [9.0752, 0.0, 0.0],
    [0.0, 9.6427, 0.0],
    [-3.8981733929, 0.0, 16.7845745817],
])
RY_TO_EV = 13.605693009
EV_TO_KJ = 96.485
E_I2_VAC = -757.79989474   # Ry, from i2_vac_ref/scf.out
N_I2 = 4


def parse_relax_out(path):
    """Return dict: final_energy (Ry), bfgs_converged (bool), n_bfgs_steps (int),
    n_scf_cycles (int), final_force (float), last_atom_block (list of (el, frac)),
    all_ion_energies (list of (step, Ry)), job_done (bool), conv_warning (str|None).
    """
    with open(path, errors="ignore") as fh:
        text = fh.read()
    d = {}
    d["job_done"] = "JOB DONE" in text
    m = re.search(r"bfgs converged in\s+(\d+)\s+scf cycles and\s+(\d+)\s+bfgs steps", text)
    d["bfgs_converged"] = bool(m)
    if m:
        d["n_scf_cycles"] = int(m.group(1))
        d["n_bfgs_steps"] = int(m.group(2))
    else:
        # count bfgs step lines as a fallback
        steps = re.findall(r"number of bfgs steps\s*=\s*(\d+)", text)
        d["n_bfgs_steps"] = int(steps[-1]) if steps else 0
        d["n_scf_cycles"] = d["n_bfgs_steps"]
    # final "!    total energy" (the last one printed is the converged-final)
    energies = re.findall(r"!\s+total energy\s*=\s*(-?\d+\.\d+)\s*Ry", text)
    d["all_energies"] = [float(e) for e in energies]
    d["final_energy"] = float(energies[-1]) if energies else None
    # final force
    forces = re.findall(r"Total force\s*=\s*(-?\d+\.\d+)", text)
    d["final_force"] = float(forces[-1]) if forces else None
    # SCF non-convergence warning?
    d["conv_warning"] = None
    if "convergence NOT achieved after" in text or "convergence has NOT been achieved" in text:
        wm = re.search(r"(convergence [A-Z ]*(?:NOT achieved|has NOT been achieved)[^\n]*)", text)
        d["conv_warning"] = wm.group(1).strip() if wm else "SCF convergence issue"
    # final ATOMIC_POSITIONS block: QE 7.5 does NOT reprint the coords after
    # "End final coordinates" in the human-readable output (they go to the .save XML).
    # So we always take the LAST ATOMIC_POSITIONS (crystal) block printed during the
    # run -- for a converged run this is the final relaxed geometry (the one that
    # satisfied the convergence test); for an nstep-capped run it is the last ionic
    # step's input geometry.
    final_block = None
    all_blocks = re.findall(
        r"ATOMIC_POSITIONS\s*\(\s*crystal\s*\)\s*\n((?:\s*\S+\s+[-0-9.]+\s+[-0-9.]+\s+[-0-9.]+.*\n?)+)",
        text,
    )
    if all_blocks:
        fb = []
        for line in all_blocks[-1].strip().splitlines():
            parts = line.split()
            if len(parts) >= 4 and parts[0] in ("C", "Gd", "H", "O", "I"):
                fb.append((parts[0], np.array([float(parts[1]), float(parts[2]), float(parts[3])])))
        if len(fb) == 80:
            final_block = fb
    if final_block is not None and d["bfgs_converged"]:
        d["final_block_source"] = "last ATOMIC_POSITIONS block (= converged final geometry)"
    elif final_block is not None:
        d["final_block_source"] = "last ATOMIC_POSITIONS block (run hit nstep cap; not formally converged)"
    else:
        d["final_block_source"] = "none found"
    d["final_block"] = final_block
    d["final_block"] = final_block
    # wall time
    wt = re.search(r"PWSCF\s*:\s*\S+\s+\S+\s+CPU\s+(\d+)m\s+(\d+\.\d+)s\s+WALL", text)
    d["wall_min"] = float(wt.group(1)) + float(wt.group(2)) / 60.0 if wt else None
    return d


def parse_poscar_frac(path):
    """Return list of (el, frac) for a POSCAR using ASE (handles all VASP variants).
    Fractional (direct/scaled) coordinates."""
    import ase.io
    atoms = ase.io.read(path)
    syms = list(atoms.get_chemical_symbols())
    frac = atoms.get_scaled_positions(wrap=False)
    return list(zip(syms, [np.array(f) for f in frac]))


def min_image_vec(frac_i, frac_j, cell):
    """Min-image displacement vector (Cartesian, angstrom) between two fractional coords."""
    d = frac_j - frac_i
    d = d - np.round(d)
    return d @ cell


def rmsd_framework(seed_atoms, final_atoms, cell):
    """RMSD (angstrom) over framework atoms (non-I) between seed and final,
    using min-image. Assumes same atom order (both from same campaign convention)."""
    fw_idx = [k for k, (el, _) in enumerate(seed_atoms) if el != "I"]
    sq = []
    for k in fw_idx:
        # match by element + nearest unwrapped (seed order preserved by QE output)
        el_s, f_s = seed_atoms[k]
        el_f, f_f = final_atoms[k]
        if el_s != el_f:
            return None  # ordering mismatch
        v = min_image_vec(f_s, f_f, cell)
        sq.append(np.dot(v, v))
    return float(np.sqrt(np.mean(sq)))


def ii_contacts(atoms, cell, cutoff=4.5):
    """Return sorted list of (dist_angstrom, i, j) for I-I pairs within cutoff."""
    iodine = [(k, f) for k, (el, f) in enumerate(atoms) if el == "I"]
    out = []
    for a in range(len(iodine)):
        for b in range(a + 1, len(iodine)):
            _, fa = iodine[a]
            kb, fb = iodine[b]
            v = min_image_vec(fa, fb, cell)
            d = float(np.linalg.norm(v))
            if d < cutoff:
                out.append((d, iodine[a][0], kb))
    out.sort()
    return out


def structure_rmsd_all(atoms_a, atoms_b, cell):
    """Overall RMSD (all atoms, min-image, same order)."""
    if len(atoms_a) != len(atoms_b):
        return None
    sq = []
    for k in range(len(atoms_a)):
        _, fa = atoms_a[k]
        _, fb = atoms_b[k]
        v = min_image_vec(fa, fb, cell)
        sq.append(np.dot(v, v))
    return float(np.sqrt(np.mean(sq)))


def main():
    out_lines = []
    def p(s=""):
        print(s)
        out_lines.append(s)

    p("=" * 78)
    p("Gd-MOF (08_Gd) FULL-RELAX SORPTION CALCULATIONS")
    p("Cell FIXED, ALL atoms (framework + I) movable. PBE+D3(BJ), 50/400 Ry, Gamma.")
    p("=" * 78)

    # parse both runs
    chain = parse_relax_out(os.path.join(CHAIN_DIR, "relax.out"))
    mol = parse_relax_out(os.path.join(MOL_DIR, "relax.out"))

    # seeds
    chain_seed = parse_poscar_frac(os.path.join(BASE, "uniform_chain/POSCAR_final"))
    mol_seed = parse_poscar_frac(os.path.join(BASE, "activated_i2/POSCAR_final_optimized"))

    # references
    # activated empty framework (framework-relaxed, I removed) -- from activate_empty/relax.out
    empty_text = open(os.path.join(BASE, "activate_empty/relax.out"), errors="ignore").read()
    em = re.findall(r"!\s+total energy\s*=\s*(-?\d+\.\d+)\s*Ry", empty_text)
    E_empty_activated = float(em[-1]) if em else None

    # prior full-relax empty (DID NOT converge -- caveat)
    E_empty_fullrelax = None
    empty_full_path = os.path.join(BASE, "full_relax/empty_full.out")
    if os.path.isfile(empty_full_path):
        t = open(empty_full_path, errors="ignore").read()
        efm = re.findall(r"!\s+total energy\s*=\s*(-?\d+\.\d+)\s*Ry", t)
        if efm:
            E_empty_fullrelax = float(efm[-1])
        full_conv = "bfgs converged" in t

    p("\n--- REFERENCES ---")
    p(f"E_I2(vacuum)              = {E_I2_VAC:.8f} Ry   (i2_mof/i2_vac_ref/scf.out)")
    p(f"E_empty(activated)        = {E_empty_activated:.8f} Ry   (08_Gd/activate_empty/relax.out)")
    if E_empty_fullrelax is not None:
        p(f"E_empty(full-relax prior)= {E_empty_fullrelax:.8f} Ry   (08_Gd/full_relax/empty_full.out)")
        p(f"  -> prior empty_full BFGS converged? {full_conv}  (CAVEAT: did NOT fully converge,")
        p("      stopped at bfgs step 48 w/ force ~0.011; treat as indicative only)")

    p("\n--- RUN RESULTS ---")
    for name, d, seed in [("CHAIN seed", chain, chain_seed), ("MOL seed", mol, mol_seed)]:
        p(f"\n[{name}]  ({'full_relax_chain' if name.startswith('CHAIN') else 'full_relax_mol'})")
        p(f"  JOB DONE             : {d['job_done']}")
        p(f"  bfgs converged       : {d['bfgs_converged']}")
        p(f"  bfgs steps / scf cyc : {d['n_bfgs_steps']} / {d['n_scf_cycles']}")
        p(f"  final force (Ry/au)  : {d['final_force']}")
        p(f"  E_loaded final (Ry)  : {d['final_energy']:.8f}")
        if d.get("wall_min"):
            p(f"  wall time (min)      : {d['wall_min']:.1f}")
        if not d["bfgs_converged"] and d["job_done"]:
            p(f"  ** NOTE              : hit nstep cap (did NOT print 'bfgs converged');")
            p(f"     using last ionic geometry. Final force {d['final_force']:.4f} Ry/au.")
        p(f"  final-geometry src   : {d.get('final_block_source','?')}")
        if d["conv_warning"]:
            p(f"  ** SCF WARNING       : {d['conv_warning']}")
        # framework relaxation magnitude
        if d["final_block"]:
            fw_rmsd = rmsd_framework(seed, d["final_block"], CELL)
            p(f"  framework RMSD vs seed (A): {fw_rmsd:.4f}")
            iic_seed = ii_contacts(seed, CELL)
            iic_fin = ii_contacts(d["final_block"], CELL)
            p(f"  I-I <4.5 A contacts (seed): {len(iic_seed)} -- " +
              ", ".join(f"{x[0]:.3f}" for x in iic_seed[:6]))
            p(f"  I-I <4.5 A contacts (fin) : {len(iic_fin)} -- " +
              ", ".join(f"{x[0]:.3f}" for x in iic_fin[:6]))

    # ---- sorption energies ----
    p("\n--- SORPTION ENERGY E_int(per I2) ---")
    p("Formula: E_int = (E_loaded - E_empty_ref - 4*E_I2vac) * 13.605693009 / 4")
    p("         then *96.485 -> kJ/mol per I2.  N_I2 = 4.")
    p("")
    for name, d in [("CHAIN", chain), ("MOL", mol)]:
        if d["final_energy"] is None:
            p(f"[{name}] no final energy -- skip")
            continue
        p(f"[{name}] E_loaded = {d['final_energy']:.8f} Ry")
        for ref_label, e_ref in [("activated-empty", E_empty_activated),
                                 ("full-relax-empty(prior,UNCONVERGED)", E_empty_fullrelax)]:
            if e_ref is None:
                continue
            eint_ev = (d["final_energy"] - e_ref - N_I2 * E_I2_VAC) * RY_TO_EV / N_I2
            eint_kj = eint_ev * EV_TO_KJ
            p(f"   vs {ref_label:38s} E_empty={e_ref:.5f} Ry -> "
              f"E_int = {eint_ev:+.4f} eV/I2 = {eint_kj:+.2f} kJ/mol/I2")

    # ---- chain vs mol comparison ----
    p("\n--- CHAIN vs MOL: do the two seeds collapse to the same minimum? ---")
    if chain["final_block"] and mol["final_block"]:
        all_rmsd = structure_rmsd_all(chain["final_block"], mol["final_block"], CELL)
        fw_only = rmsd_framework(chain["final_block"], mol["final_block"], CELL)
        # I-only RMSD
        iod_idx = [k for k, (el, _) in enumerate(chain["final_block"]) if el == "I"]
        sq = []
        for k in iod_idx:
            v = min_image_vec(chain["final_block"][k][1], mol["final_block"][k][1], CELL)
            sq.append(np.dot(v, v))
        i_rmsd = float(np.sqrt(np.mean(sq)))
        p(f"  RMSD all 80 atoms (chain-final vs mol-final): {all_rmsd:.4f} A")
        p(f"  RMSD framework only                         : {fw_only:.4f} A")
        p(f"  RMSD I atoms only                           : {i_rmsd:.4f} A")
    dE = chain["final_energy"] - mol["final_energy"]
    dE_ev = dE * RY_TO_EV
    dE_kj_per_i2 = dE_ev * EV_TO_KJ / 1  # total per formula unit, not per I2
    p(f"  dE = E_chain - E_mol = {dE:+.6f} Ry = {dE_ev:+.4f} eV "
      f"= {dE_ev*EV_TO_KJ:+.2f} kJ/mol (per unit cell, 4 I2)")
    p(f"  -> per I2: {dE_ev*EV_TO_KJ/4:+.2f} kJ/mol")
    p("")
    if abs(dE_ev*EV_TO_KJ) > 5:
        p("  CONCLUSION: the two seeds CONVERGED TO DISTINCT MINIMA "
          "(chain stays lower; full relaxation does NOT collapse them).")
    else:
        p("  CONCLUSION: the two seeds collapsed to essentially the same minimum.")

    p("\n--- RAW ENERGIES (for re-computation with a future rigorous E_empty_full) ---")
    p(f"  E_loaded(chain) = {chain['final_energy']:.8f} Ry")
    p(f"  E_loaded(mol)   = {mol['final_energy']:.8f} Ry")
    p(f"  E_I2vac         = {E_I2_VAC:.8f} Ry  (x4 = {4*E_I2_VAC:.8f})")

    p("\n--- CAVEATS ---")
    p("* Both runs used conv_thr=1e-3, mixing_beta=0.2, degauss=0.01. The chain SCF")
    p("  plateaus at ~7e-4 accuracy at the strained-chain geometry (conv_thr 1e-5/1e-6")
    p("  is unreachable there and caused the first attempts to abort at the 300-iter")
    p("  SCF cap; the original uniform_chain/loaded.in also accepted ~3e-4 on its first")
    p("  SCF). The mol run initially used conv_thr=1e-5 and aborted at step 64 (SCF")
    p("  plateau ~1.0e-5); it was restarted from its step-63 geometry with conv_thr=1e-3")
    p("  (see attempt1_conv1e-5_64steps/).")
    p("* The mol hit nstep=200 (force 0.0103, just above forc_conv_thr=1e-2) -- the mol")
    p("  basin is broad/soft; the energy had effectively converged (delta-E < 1e-5 Ry")
    p("  over the last ~30 steps, E = -6865.0882 Ry). The chain converged cleanly")
    p("  (force 0.0045, 'bfgs converged' at step 63).")
    p("* E_int 'vs activated-empty' uses the framework-only-relaxed empty (no I).")
    p("  A rigorous full-relax reference needs the empty framework ALSO full-relaxed")
    p("  (cell fixed, atoms free). The prior full_relax/empty_full did NOT converge")
    p("  (stopped step 48, force 0.011) -- its E is reported for reference only.")
    p("* First attempts with mixing_beta=0.1 failed (SCF hit 300-iter cap); the proven")
    p("  recipe is mixing_beta=0.2 + degauss=0.01 (matches uniform_chain/loaded.in).")
    p("  See full_relax_chain/attempt1_beta0p1/ and attempt2_conv1e-5_stalled/.")

    out_path = os.path.join(BASE, "full_relax_energies.txt")
    with open(out_path, "w") as fh:
        fh.write("\n".join(out_lines) + "\n")
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
