#!/usr/bin/env python3
"""Compare the mixed (POSCAR_final) packing vs the two uniform variants for the 3 pts systems.

  E_int    = E_loaded - E_empty - N_I2 * E(I2_vac)     [eV]
  E_per_I2 = E_int / N_I2

E_empty is reused from each parent's empty_mof.out (same framework, frozen).
E_mixed baseline = the lowest of each parent's loaded relaxes (loaded_seeded for 13_Y,
loaded_tight for 06_Tb/08_Gd).
"""
import os, glob
RY = 13.605693009
E_I2 = -757.80004425
BASE = os.path.dirname(os.path.abspath(__file__))
STRUCTS = {"13_Y": 4, "06_Tb": 4, "08_Gd": 4}
MIXED = {"13_Y": "loaded_seeded", "06_Tb": "loaded_tight", "08_Gd": "loaded_tight"}


def e(path):
    if not os.path.isfile(path):
        return None, False
    val = None
    txt = open(path, errors="ignore").read()
    for line in txt.splitlines():
        if "!    total energy" in line:
            val = float(line.split("=", 1)[1].split()[0])
    return val, ("JOB DONE" in txt)


def main():
    print(f"{'system':7} {'variant':17} {'E_loaded(Ry)':>15} {'E_int(eV)':>11} "
          f"{'E/I2(eV)':>10} {'E/I2(kJ/mol)':>13}  status")
    for s in sorted(STRUCTS):
        n = STRUCTS[s]
        eE, _ = e(f"{BASE}/{s}/empty_mof.out")
        rows = []
        for label, sub in [("mixed", MIXED[s]), ("uniform_chain", "uniform_chain/loaded"),
                           ("uniform_herring", "uniform_herring/loaded")]:
            p = f"{BASE}/{s}/{sub}.out" if sub.endswith("loaded") else f"{BASE}/{s}/{sub}.out"
            el, dl = e(p)
            if eE is not None and el is not None:
                ei = (el - eE - n * E_I2) * RY
                rows.append((s, label, el, ei, ei / n, ei / n * 96.485, dl))
            else:
                rows.append((s, label, el, None, None, None, dl))
        # marker: best (lowest E/I2)
        best = min((r for r in rows if r[4] is not None), default=None, key=lambda r: r[4])
        for r in rows:
            mark = "  <== lowest" if (best and r[1] == best[1]) else ""
            st = "OK" if r[6] else ("running" if r[2] else "no-file")
            print(f"{r[0]:7} {r[1]:17} {(r[2] or 0):>15.6f} {(r[3] or 0):>11.4f} "
                  f"{(r[4] or 0):>10.4f} {(r[5] or 0):>13.2f}  {st}{mark}")
        print()
    print(f"E(I2_vac) = {E_I2} Ry   E_empty reused from parent empty_mof.out")


if __name__ == "__main__":
    main()
