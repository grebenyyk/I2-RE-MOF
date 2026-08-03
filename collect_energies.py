#!/usr/bin/env python3
"""
Collect total energies from QE outputs and compute I2-framework interaction energies.

  E_int    = E(loaded) - E(empty) - N_I2 * E(I2_vac)
  E_per_I2 = E_int / N_I2   (eV)

E(I2_vac) is the shared isolated-I2 reference from i2_ce-uio-66/01_i2_vac
(PBE+D3(BJ), 50/400, I.pbe-n-kjpaw_psl): -757.80004425 Ry.

Dispersion contribution: from the "DFT-D3 Dispersion =" term in each energy breakdown,
  E_disp_int = D(loaded) - D(empty) - N_I2 * D(I2_vac)   (Ry)
reported per I2 in kJ/mol. The remainder (E_int - E_disp_int) is the DFT (electronic)
part of the binding.

Reads <system>/empty_mof.out and <system>/loaded.out (prefers loaded_tight.out if present).
A system is "done" only if its output contains JOB DONE; otherwise the energy shown is the
latest SCF value (interim) and E_int is flagged as not-final.
"""
import os, glob

RY_TO_EV = 13.605693009
EV_TO_KJMOL = 96.48533212  # 1 eV per molecule = 96.485 kJ/mol
E_I2_VAC = -757.80004425  # Ry
I2_VAC_OUT = "/workspace/dimitry/calc/qe/i2_ce-uio-66/01_i2_vac/i2_vac.out"

# system label -> N_I2 (I atoms / 2), from the POSCARs
N_I2 = {"13_Y": 4, "08_Gd": 4,
        "05_Cd": 5, "07_Cd": 6, "10_Cd": 5, "17_Cd": 2,
        "06_Tb": 4, "12_Tb": 6, "16_U": 11, "18_U": 10}


def _last_match(path, needle):
    if not os.path.isfile(path):
        return None
    val = None
    for line in open(path, errors="ignore"):
        if needle in line and "=" in line:
            try:
                val = float(line.split("=", 1)[1].split()[0])
            except (ValueError, IndexError):
                pass
    return val


def final_energy(outpath):
    """Return (energy_Ry, done_bool). Energy = last '! total energy' line."""
    e = _last_match(outpath, "!    total energy")
    done = os.path.isfile(outpath) and "JOB DONE" in open(outpath, errors="ignore").read()
    return e, done


def disp_energy(outpath):
    """Last 'DFT-D3 Dispersion =' value in Ry (None if absent)."""
    return _last_match(outpath, "DFT-D3 Dispersion")


def best_loaded(d):
    """Prefer loaded_seeded.out > loaded_tight.out > loaded.out. Returns (e, done, path)."""
    for name in ("loaded_seeded.out", "loaded_tight.out", "loaded.out"):
        p = os.path.join(d, name)
        e, done = final_energy(p)
        if e is not None:
            return e, done, p
    return None, False, None


def main():
    base = os.path.dirname(os.path.abspath(__file__))
    d_i2vac = disp_energy(I2_VAC_OUT)
    rows = []
    for d in sorted(glob.glob(os.path.join(base, "[01]?_*"))):
        label = os.path.basename(d)
        if label not in N_I2:
            continue
        n = N_I2[label]
        eE, dE = final_energy(os.path.join(d, "empty_mof.out"))
        eL, dL, pL = best_loaded(d)
        dE_disp = disp_energy(os.path.join(d, "empty_mof.out"))
        dL_disp = disp_energy(pL) if pL else None
        if eE is not None and eL is not None:
            e_int = (eL - eE - n * E_I2_VAC) * RY_TO_EV
            e_per = e_int / n
        else:
            e_int = e_per = None
        if dE_disp is not None and dL_disp is not None and d_i2vac is not None:
            disp_int_kj = (dL_disp - dE_disp - n * d_i2vac) * RY_TO_EV * EV_TO_KJMOL
            disp_per_kj = disp_int_kj / n
        else:
            disp_per_kj = None
        rows.append((label, n, eE, dE, eL, dL, e_int, e_per, disp_per_kj, pL))

    print(f"{'system':8s} {'N_I2':>4s} {'E_int(eV)':>10s} {'E_int(kJ/mol)':>13s} "
          f"{'E/I2(eV)':>9s} {'E/I2(kJ/mol)':>13s} {'E_disp/I2(kJ/mol)':>17s}  status")
    for label, n, eE, dE, eL, dL, e_int, e_per, disp_per_kj, pL in rows:
        st = ("OK" if (dE and dL) else ("empty done, loaded running" if dE else "running"))
        st += f"  [{os.path.basename(pL) if pL else '?'}]"
        ei_kj = (e_int * EV_TO_KJMOL) if e_int is not None else 0
        ep_kj = (e_per * EV_TO_KJMOL) if e_per is not None else 0
        dk = f"{disp_per_kj:>17.2f}" if disp_per_kj is not None else f"{'--':>17s}"
        print(f"{label:8s} {n:>4d} {(e_int or 0):>10.4f} {ei_kj:>13.2f} "
              f"{(e_per or 0):>9.4f} {ep_kj:>13.2f} {dk}  {st}")
    print(f"\nE(I2_vac) = {E_I2_VAC} Ry ; D(I2_vac) = {d_i2vac} Ry  (shared reference)")


if __name__ == "__main__":
    main()
