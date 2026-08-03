# QTAIM plane maps from QE cubes (Multiwfn)

Publication figures: electron-density contours + gradient-line bond paths +
auto-marked critical points, rendered by **Multiwfn** from the existing
`val.cube` (valence density) files. No cubes were regenerated.

## Files
| file | what |
|------|------|
| `Y_act_gradient.png` / `Y_act_panelA.pdf` | **Panel A** — Y_act covalent I₂ (plane I76–I77–C6; I76–I77 = 2.692 Å) |
| `Y_deep_gradient.png` / `Y_deep_panelB.pdf` | **Panel B** — Y_deep polyiodide (plane I75–I74–I76; ~2.93 Å, central I74) |
| `QTAIM_two_panel.png` / `QTAIM_two_panel.pdf` | two-panel composite (a)/(b) |

Each map shows: black density contour lines; gradient-line trajectories = bond
paths; atoms labelled (element symbol); blue sticks = empirical bonds; **critical
points marked** by Multiwfn (brown = nuclear, blue = BCP (3,−1), orange = ring
(3,+1), green = cage (3,+3)); axes in Å.

## How to regenerate / tweak
```bash
cd /workspace/dimitry/calc/qe/i2_mof/qtaim/campaign/bcp_plots/multiwfn_out
# settings.ini (in THIS folder) must contain: isilent=1, iuserfunc=-3, graphformat=png
BIN=/workspace/dimitry/calc/qe/i2_mof/tools/multiwfn/Multiwfn_2026.7.11_bin_Linux/Multiwfn
export MultiwfnPATH="$(dirname $BIN)" OMP_STACKSIZE=128m
# args: cube, 3 atoms defining the plane, output name
CUBE=.../Y_act/val.cube; ATOMS=76,77,6
printf '%s\n\n4\n100\n6\n300,300\n4\n%s\n8\n3\n-8\n0\n-5\nq\n' "$CUBE" "$ATOMS" \
  | timeout 200 xvfb-run -a "$BIN"   # writes dislin.png; cp to <name>.png
```
Change map type `6`→`2` for pure contour (no gradient lines), `1` for colour-filled.
Change atom triple to reorient the plane. Flip `graphformat` to `pdf` for vector masters.

## Plane / atom choice rationale
Iodines are atoms **73–80** (8 per cell). From cube geometry:
- Y_act: I76–I77 = 2.692 Å (the only true covalent I₂ besides I79–I80); plus a
  polyiodide cluster {I73,I74,I75,I78} (I74–I75 = 2.951 Å, I73–I74 = 3.017 Å).
- Y_deep: no isolated I₂ (all I–I 2.9–3.4 Å); two {I4}-like units around I74
  (bonds to I73/I75/I76) and I78 (bonds to I77/I79/I80). Panel shows I75–I74–I76.

See also: `compose_two_panel.py` (composite builder) and project memory
`i2-mof-multiwfn-qtaim-plot-recipe`.

## Getting CP DOTS (the cube limitation + the fix)

Multiwfn **cannot locate critical points on a cube** — its grid CP-search returns
0 CPs (confirmed for the I76–I77 covalent bond). The gradient map draws
bond-path *lines* but no CP *dots*. Dots only appear once CPs are in memory.

**Fix = load critic2's already-found CPs** (`bcp.out` in each `Y_*/` dir):
1. Convert: `python3 bcpout_to_cpstxt.py <bcp.out> <val.cube> CPs.txt [--with-ring-cage]`
   (fractional→Bohr Cartesian; skips nuclear CPs; drops <0.55 Å near-core artifacts).
2. In Multiwfn: load cube → `2` (topology) → `-4` → `5` (load `CPs.txt`, ENTER) →
   `0` → `-10` (main) → `4` (plane) → `100` → `6` → grid → `4` → atoms.
3. **Telltale it worked**: post-menu now lists `4 Set details of plotting critical
   points and paths` and `6 Generate and show interbasin paths` (absent without CPs).
4. Post-menu `4` adjusts CP dot size/color; `6` draws true interbasin bond paths.

`CPs.txt` format (1 line count, then `idx x y z type`, Bohr, type 1=NCP 2=BCP 3=ring 4=cage).
Verified: loaded CP 284 (I74–I75 BCP) plots exactly between cube atoms 74,75 at 2.95 Å.
