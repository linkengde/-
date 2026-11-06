# P1I final report: 3× exploratory arc pressure and Marangoni run

## Outcome

The run formed and then solidified an Ag pool, but **did not complete the requested 3 ps hold**. The 2100 K trigger was reached at 0.62344771 ps (measured skeleton maximum 2117.8 K). At 0.80095092 ps, after 0.17750 ps of the hold, the predictive property guard shut off heat and pressure because the next off-ramp step was predicted to reach 2300.32 K. The measured skeleton maximum was 2292.3 K, below the 2300 K table ceiling. Cooling then reached the solver's no-liquid solidification criterion at 1.10183585 ps.

The peak-pool surface-height diagnostic remained slightly raised relative to the initial surface: `source_center_depression_nm = -0.2059 nm` at heat-off and `-0.2185 nm` after solidification. In this diagnostic, negative means the surface is above its initial height. The run therefore produced **no resolved permanent central pit**. The final offset is smaller than one cell (Δx = 0.30186 nm), so it is not a resolved crater-scale shape change.

## Model and run settings

- Ag70 geometry, fixed grid 64 × 64 × 76, open ambient-pressure top, periodic x/y, Δx = 0.301857485 nm; sample geometry unchanged.
- Existing single-temperature enthalpy, phase-change and composite heating implementation; reference source power 8.931899 mW.
- Marangoni enabled at 3× for tangential thermal acceleration; normal capillary stress unchanged.
- Periodic Gaussian arc-pressure profile, 501 MPa configured peak, 5 nm radius, 0.02 ps ramp. Maximum effective interface pressure recorded was 500.416 MPa; integrated z force at heat-off was −2.4573 × 10⁻⁸ N.
- Ti₃SiC₂ property table 300–2300 K. The predictive 2300 K protection, hard ceiling and divergence guard remained active.
- Before the 2100 K trigger, the source could also stop after the input-to-bottom-sink difference stayed within 0.01 mW for 10 consecutive steps. That condition did not occur; the difference was about 6.01 mW at heat-off.
- Cooling stopped when liquid Ag was absent for two consecutive steps. The final solidified state was still warm (skeleton maximum 1109.0 K), not cooled to 300 K.

All pressure values and the threefold Marangoni multiplier are exploratory model inputs, not calibrated nanoscale arc conditions or a quantitative prediction of a physical contact.

## Peak pool and conservation

At 0.80095092 ps:

| Quantity | Result |
|---|---:|
| Total liquid Ag volume | 190.949 nm³ |
| Largest connected liquid cluster | 175.242 nm³ (6,403 cells) |
| Largest-cluster bounding size | 11.772 × 11.772 × 3.924 nm |
| Skeleton / Ag maximum temperature | 2292.3 K / 2265.7 K |
| Maximum velocity at this history sample | 144.29 m/s |
| Center surface-height diagnostic | −0.2059 nm (raised) |

At final solidification, liquid Ag volume was zero, skeleton and Ag maxima were 1109.0 K and 1104.2 K, respectively, and the center surface-height diagnostic was −0.2185 nm (still raised). Maximum absolute Ag-volume drift over the run was 0.00437%. The largest absolute energy residual was 9.91 × 10⁻²² J.

## Archived outputs

`results/full_3x_hold3ps_20261106/` contains the full time-series CSVs, pressure history, pool dimensions, run configuration, logs, code hashes, initial/hold-start/heat-off/solidified field and pressure frames, final solidified dump, rendered 3D views and SHA-256 manifest. The peak-pool image shows solid Ag in light blue, Ti₃SiC₂ in dark blue, partially liquid Ag in yellow, and liquid-dominant Ag in orange.

The final dump and field snapshots were validated before the remaining periodic restart checkpoints were removed to recover workspace storage. The archived run is therefore a result archive, not a restartable run directory.

The 0.16 ps preliminary run and setup details are documented separately in `P1I_PRELIMINARY_PREFLIGHT_REPORT.md`. This full run supersedes it for the 3× parameter case.

## Storage cleanup

After the run was archived and pushed, 51 restart-checkpoint and stale-pointer files from previously completed P1G smoke runs were removed, freeing 0.954 GiB. Final dumps, logs, CSV summaries, and field snapshots were retained. Failed, ambiguous, and active runs were left untouched. The removed-file list and pre-removal SHA-256 values are recorded in `results/p1g_completed_smoke_checkpoint_cleanup_20261106.csv`.
