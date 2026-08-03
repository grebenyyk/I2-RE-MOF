# QTAIM + charge-transfer study — I2@MOF (PBE)

Nine structures: solvent-templated (st), activated (act), variable-cell (vc) x Y, Tb, Gd.

## cubes/  (9 systems x 6 cubes)
  {sys}_rho.cube       = rho17.cube — plot_num=17, PAW-reconstructed valence (less
                         noisy) -> used for BCP topology (critic2)
  {sys}_val.cube       = val.cube — plot_num=0, raw valence (loaded) -> constituent of deltarho
  {sys}_empty.cube     = empty_val.cube — empty framework density
  {sys}_i2only.cube    = i2only_val.cube — isolated I2 density
  {sys}_deltarho.cube  = charge-density difference = val - empty - i2only  (verified exactly)
  {sys}_elf.cube       = ELF (plot_num=8)
  rho (plot_num=17) and val (plot_num=0) are the same physical density, two pp.x
  reconstructions; plot_num=17 preferred for topology (PAW reconstruction smooths
  the valence field), plot_num=0 used to build the charge-density difference.

## critic2_outputs/
  {sys}_bcp_raw.out        = critic2 `auto seed pair dist 4.0` output — the CANONICAL
                             CP search; locates I-I, I-O AND Y/Ln..I BCPs. The default
                             pair-seeding distance misses the longer Y/Ln..I contacts
                             (3.2-3.4 A); raising the cutoff to 4.0 A captures them.
  {sys}_bcp_filtered.out   = (3,-1) BCPs after geometry filter (d1<2.6 A, d2/d1<2)
  {sys}_bcp.cri            = critic2 input (original default-pair run; superseded)
  {sys}_bader_charges.out  = ct.out — Bader/integrated charges (charge-transfer table)

## scripts/  — seed_bcp.py (cube Newton/interp), bcp_energy.py (Abramov G/V/H),
   make_dedup_ii/io.py (I-I, I-O table builders), make_mi_native.py (M-I table
   from the seed-pair-4.0 critic2 output), profile_bcp.py, bcpout_to_cpstxt.py.

## tables/  — BCP_II_full (75 I-I, both chains), BCP_IO_full (69 I-O),
   BCP_MI_native (10 Y/Ln..I).  .docx + .csv each.

## System -> source dir (i2_mof/qtaim/campaign/)
Y-MOFst=Y_str  Tb-MOFst=Tb_str  Gd-MOFst=Gd_str
Y-MOFact=Y_actc Tb-MOFact=Tb_act Gd-MOFact=Gd_act
Y-MOFvc=Y_deep  Tb-MOFvc=Tb_deep Gd-MOFvc=Gd_deep

## Methods recap
BCPs located with critic2 (v1.2.1005) on the PAW-reconstructed valence cube
(plot_num=17) using pair-seeding to 4.0 A, so that the longer Y/Ln..I contacts
(3.2-3.4 A) — which exceed the default pair-seeding distance — are captured
natively. Spurious near-core (3,-1) discarded (d1<2.6 A, d2/d1<2). I-I
rho>=0.02 e/bohr3; I-O (side-on, ~90 deg) retained <4 A; Y/Ln..I ionic
(rho~0.02 st, ~0.013 vc). G/V/H: G=(3/10)(3pi^2)^(2/3)rho^(5/3)+(1/6)nabla2rho;
V=(1/4)nabla2rho-2G; H=G+V (a.u.). Charge transfer: Bader integration of
deltarho = val - empty - i2only (plot_num=0).
