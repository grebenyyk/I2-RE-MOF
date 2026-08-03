# Raman band → I–I contact assignment (I sublattice, PBE+D3)

Frequencies from the finite-difference Hessian (framework frozen, Δ=0.01 Å). Relative activities and contact percentages from the bond-polarizability model (dα/dr ∝ r⁻⁴, r₀=2.67 Å gas I₂), applied only to chemically bonded I–I pairs (≤3.25 Å — the natural coordination gap; see geometry table). Non-bonded van-der-Waals contacts (≥3.25 Å, e.g. the 3.61 Å branching-I···in-chain distance) carry no bond polarizability and are excluded: they would otherwise leak spurious intensity into chain-bending modes. Mode→contact fractions are read off the saved eigenvectors; the spectrum and this table share the same eigenvectors/model.


## I–I contact geometry (< 4.0 Å, min-image). ≥3.25 Å = non-bonded vdw (excluded from Raman model)

| variant | contact | r (Å) | motif |
|---|---|---|---|
| Y-MOFst | I03-I04 | 2.93 | chain bond |
| Y-MOFst | I07-I08 | 2.93 | chain bond |
| Y-MOFst | I06-I07 | 2.95 | chain bond |
| Y-MOFst | I02-I03 | 2.95 | chain bond |
| Y-MOFst | I01-I02 | 2.99 | chain bond |
| Y-MOFst | I05-I06 | 2.99 | chain bond |
| Y-MOFst | I02-I04 | 3.12 | secondary bond |
| Y-MOFst | I06-I08 | 3.12 | secondary bond |
| Y-MOFst | I01-I04 | 3.61 | vdw contact (non-bonded) |
| Y-MOFst | I05-I08 | 3.61 | vdw contact (non-bonded) |
| Y-MOFact | I05-I06 | 2.74 | I2-like pair |
| Y-MOFact | I01-I02 | 2.74 | I2-like pair |
| Y-MOFact | I03-I04 | 2.89 | chain bond |
| Y-MOFact | I07-I08 | 2.89 | chain bond |
| Y-MOFact | I06-I07 | 3.01 | secondary bond |
| Y-MOFact | I02-I04 | 3.02 | secondary bond |
| Y-MOFact | I06-I08 | 3.03 | secondary bond |
| Y-MOFact | I02-I03 | 3.04 | secondary bond |
| Y-MOFvc | I01-I05 | 2.82 | chain bond |
| Y-MOFvc | I02-I03 | 2.93 | chain bond |
| Y-MOFvc | I06-I07 | 2.94 | chain bond |
| Y-MOFvc | I06-I08 | 2.94 | chain bond |
| Y-MOFvc | I02-I04 | 2.94 | chain bond |
| Y-MOFvc | I07-I08 | 2.94 | chain bond |
| Y-MOFvc | I03-I04 | 2.94 | chain bond |
| Y-MOFvc | I05-I06 | 3.40 | vdw contact (non-bonded) |
| Y-MOFvc | I01-I02 | 3.42 | vdw contact (non-bonded) |

## Band assignment (spectral peaks; broadened intensity + motif)

`rel. int.` = broadened peak height as a % of each variant's tallest peak (i.e. exactly the height in Fig. `raman_Y_st_act_vc`). `n` = number of near-degenerate modes in the peak — a 4-mode cluster towers over a single mode of equal per-mode activity.

| variant | ν (cm⁻¹) | n | rel. int. | dominant motif + top contacts |
|---|---|---|---|---|
| Y-MOFst | 152--156 | 4 | 100 | chain bond (84%): I02-I03 (2.95 Å; 17%), I06-I07 (2.95 Å; 17%), I07-I08 (2.93 Å; 17%) |
| Y-MOFst | 126 | 2 | 37 | chain bond (79%): I01-I02 (2.99 Å; 31%), I05-I06 (2.99 Å; 31%), I06-I08 (3.12 Å; 11%) |
| Y-MOFact | 160 | 2 | 86 | chain bond (58%): I03-I04 (2.89 Å; 31%), I07-I08 (2.89 Å; 27%), I05-I06 (2.74 Å; 13%) |
| Y-MOFact | 150--152 | 2 | 100 | I2-like pair (82%): I01-I02 (2.74 Å; 44%), I05-I06 (2.74 Å; 39%), I07-I08 (2.89 Å; 7%) |
| Y-MOFact | 144--145 | 2 | 80 | secondary bond (96%): I06-I08 (3.03 Å; 29%), I02-I03 (3.04 Å; 27%), I02-I04 (3.02 Å; 21%) |
| Y-MOFvc | 173--176 | 4 | 100 | chain bond (100%): I02-I03 (2.93 Å; 17%), I06-I08 (2.94 Å; 17%), I06-I07 (2.94 Å; 17%) |
| Y-MOFvc | 159 | 1 | 33 | chain bond (100%): I01-I05 (2.82 Å; 99%) |
