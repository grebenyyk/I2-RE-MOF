# Landscape 3-variants archive

Configs (HSE06-D3/4x4x2 single-points on PBE/4x4x2-relaxed geometries):

| Variant | Motif | Config | E (Ry) | dEads (kJ/mol/I2) | Lowest |
|---|---|---|---|---|---|
| st | CIF | Y-MOFst_polyiodide-CIF | -5766.47005 | -130.08 |  |
| st | poly | Y-MOFst_polyiodide-S1-clustered | -5766.52530 | -148.21 | LOWEST polyiodide (st) |
| st | poly | Y-MOFst_polyiodide-S2 | -5766.49377 | -137.86 |  |
| st | poly | Y-MOFst_polyiodide-S3 | -5766.50425 | -141.30 |  |
| st | poly | Y-MOFst_polyiodide-S4 | -5766.48588 | -135.27 |  |
| st | poly | Y-MOFst_polyiodide-S5 | -5766.49968 | -139.80 |  |
| st | poly | Y-MOFst_polyiodide-S6 | -5766.49422 | -138.01 |  |
| st | mol | Y-MOFst_molecular-S2 | -5766.34319 | -88.45 |  |
| st | mol | Y-MOFst_molecular-S3 | -5766.34320 | -88.45 |  |
| st | mol | Y-MOFst_molecular-S4 | -5766.33780 | -86.67 |  |
| st | mol | Y-MOFst_molecular-S6 | -5766.36941 | -97.05 | LOWEST molecular (st) |
| act | CIF | Y-MOFact_polyiodide-CIF | -5767.20025 | -69.64 |  |
| act | poly | Y-MOFact_polyiodide-S1 | -5767.21388 | -74.11 |  |
| act | poly | Y-MOFact_polyiodide-S2 | -5767.21640 | -74.93 | LOWEST polyiodide (act) |
| act | poly | Y-MOFact_polyiodide-S3 | -5767.21227 | -73.58 |  |
| act | poly | Y-MOFact_polyiodide-S4 | -5767.20474 | -71.11 |  |
| act | poly | Y-MOFact_polyiodide-S5 | -5767.21590 | -74.77 |  |
| act | poly | Y-MOFact_polyiodide-S6 | -5767.21512 | -74.51 |  |
| act | mol | Y-MOFact_molecular-S1 | -5767.25107 | -86.31 | LOWEST molecular (act) |
| act | mol | Y-MOFact_molecular-S2 | -5767.16673 | -58.63 |  |
| act | mol | Y-MOFact_molecular-S3 | -5767.16672 | -58.63 |  |
| act | mol | Y-MOFact_molecular-S4 | -5767.23574 | -81.28 |  |
| act | mol | Y-MOFact_molecular-S6 | -5767.01046 | -7.35 |  |
| vc | CIF | Y-MOFvc_polyiodide-CIF | -5767.32487 | -111.98 |  |
| vc | poly | Y-MOFvc_polyiodide-S1-clustered | -5767.34440 | -118.39 | LOWEST polyiodide (vc) |
| vc | poly | Y-MOFvc_polyiodide-S2 | -5767.32805 | -113.03 |  |
| vc | poly | Y-MOFvc_polyiodide-S3 | -5767.33459 | -115.17 |  |
| vc | poly | Y-MOFvc_polyiodide-S4 | -5767.31952 | -110.23 |  |
| vc | poly | Y-MOFvc_polyiodide-S5 | -5767.31668 | -109.30 |  |
| vc | poly | Y-MOFvc_polyiodide-S6 | -5767.33147 | -114.15 |  |
| vc | mol | Y-MOFvc_molecular-S1 | -5767.14877 | -54.19 |  |
| vc | mol | Y-MOFvc_molecular-S2 | -5766.97695 | +2.20 |  |
| vc | mol | Y-MOFvc_molecular-S3 | -5766.98725 | -1.18 |  |
| vc | mol | Y-MOFvc_molecular-S4 | -5767.15157 | -55.11 | LOWEST molecular (vc) |
| vc | mol | Y-MOFvc_molecular-S5 | -5767.11624 | -43.51 |  |

The 6 lowest (min polyiodide + min molecular per variant) are also copied as `<name>__LOWEST.cif` for quick identification.

References: rd_gas, rd_empty_Y{st,act,vc} (gas I2 + empty frameworks).
Gas = -757.80001 Ry.
