# r2SCAN cross-check — ΔEads (TOTAL, includes D3), kJ/I2
CORRECTED 2026-08-07: prior version double-counted D3 (bug: le+ld where !E already incl D3).
QE convention (electrons.f90: etot=etot+edftd3): "! total energy" INCLUDES D3.
=> total = (!E_load - !E_empty - 4 !E_gas); D3 = (D3c_load - D3c_empty - 4 D3c_gas); elec = total - D3.

## Total ΔEads (kJ/I2)
| system | PBE-PAW prod (e50/400) | PBE-NC (e70/560) | r2SCAN-NC relaxed (e70/560) |
|---|---|---|---|
| Yact | -67.1  | -57.9  | -53.3  |
| Yst  | -136.5 | -138.9 | -139.9 |

=> r2SCAN-NC AGREES with PBE: Yst within 3 kJ, Yact within 14 kJ (9 ecut e50->e70 + 5 functional).
   The cross-check VALIDATES the production PBE-PAW binding energies. NO 56 kJ gap (that was the double-count bug).

## electronic / D3 decomposition (functionals agree on total, redistribute)
| system-func | total | D3 | elec |
|---|---|---|---|
| Yact PBE-PAW | -67.1 | -64.9 | -2.2 |
| Yact r2SCAN-NC | -53.3 | ~-35 | ~-18 |
| Yst PBE-PAW | -136.5 | -62.4 | -74.1 |
| Yst r2SCAN-NC | -139.9 | ~-35 | ~-105 |
(r2SCAN: more electronic, less D3 than PBE — D3 params functional-specific. Same total.)

## PP-cross (valid; D3 cancels): PAW-vs-NC <0.5 kJ for I AND framework.
## ecut: e70 vs e90 (NC, Yact) = +0.28 kJ; Yst e50 vs e70 ~2 kJ; Yact e50 vs e70 ~9 kJ (empty-shift, minor).
