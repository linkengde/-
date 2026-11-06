# P1I 6 mW composite-source run

## Outcome

The Ag70 run used the existing composite Ag–Ti₃SiC₂ heat-source distribution, renormalized from its 8.931899 mW reference integral to **6.000 mW**. It started at 300 K, reached the 2100 K skeleton target at 2.276773 ps, then continued heating with 3× tangential Marangoni forcing and a configured 501 MPa Gaussian arc-pressure load.

The requested 3 ps post-target hold did not complete. At restart from a valid same-run checkpoint at 3.294562 ps, the pre-step skeleton peak was 2276.01 K, already above the new 2250 K pre-step shutdown threshold. The solver therefore shut off the heat source and pressure before the next thermal solve, at 3.296380 ps, after **1.019608 ps** of the hold. The hard 2300 K property-table ceiling was not crossed in this accepted trajectory. Cooling continued to the solver's two-step no-liquid criterion at 3.878626 ps.

An earlier attempt using only the predictive guard stopped at 3.514086 ps after the skeleton briefly reached 2300.83177 K. That out-of-table state was rejected as a restart source; its log and dump hash are retained under `results/full_6mW_3x_safe_20261106/diagnostics/`. The accepted continuation starts from the complete 3.294562 ps checkpoint, before the failed step.

## Pool result

The largest recorded total liquid-Ag volume was **475.187 nm³ at 3.304821 ps**, shortly after the heat/pressure cutoff ramp began. The largest connected majority-liquid cluster at the same saved statistic time was **444.568 nm³** across 16,337 cells, with the reported bounding size **17.508 × 17.206 × 5.433 nm**. These are continuum-grid diagnostics on the 64 × 64 × 76 grid (Δx = 0.301857 nm).

At solidification, liquid Ag was zero to the solver tolerance. The maximum accepted time-series temperature was 2274.55 K; the separate pre-step guard read 2276.01 K, while the heat-off saved field was 2274.27 K. The maximum effective interface pressure recorded was 500.416 MPa. The final Ag-volume drift was −0.06995%; the largest absolute energy residual was 3.46 × 10⁻²¹ J (2.43 × 10⁻⁵% of cumulative input energy at its worst relative point). The final sampled surface-height diagnostic remained slightly raised relative to its initial reference (mean −0.2092 nm; source-center −0.1646 nm, where negative denotes upward displacement), so this run did not form a resolved permanent central pit.

## Configuration and limits

- Existing Ag70 geometry; no geometry or material-mask changes.
- Single-temperature enthalpy model with the existing phase change and bottom heat sink.
- Composite Ag + Ti₃SiC₂ source shape preserved and normalized to 6.000 mW.
- 3× tangential Marangoni multiplier; normal capillary stress unchanged.
- 501 MPa configured peak periodic-Gaussian arc pressure, 5 nm radius; pressure was switched off with the heat source.
- Ti₃SiC₂ table range 300–2300 K. The accepted saved-field maximum was 2274.27 K; the pre-step shutdown check read 2276.01 K. No temperature clipping was used.
- Ag volume conservation and energy residuals are reported in the merged full-run CSVs.

## Files

The result archive is `results/full_6mW_3x_safe_20261106/`.

- `p1i_6mW_peak_pool_3d.png`: three-dimensional phase view at heat-off.
- `p1i_6mW_final_solidified_3d.png`: final solidified geometry.
- `p1i_6mW_temperature_pool_history.png`: temperature, liquid volume, and power history.
- `p1h_history.csv`: merged 300 K start through solidification time series; the restart join is preserved in the manifest.
- `p1h_peak_pool_dimensions.csv`: largest connected liquid-cluster statistics.
- `p1h_arc_pressure_history.csv`: applied source scaling and interface pressure history.
- `fields/`: initial, target, hold, heat-off, and solidified temperature/material/velocity/pressure snapshots.
- `run_manifest.json`, `run_config.txt`, `code_hashes_*.sha256`, and `SHA256SUMS`: reproducibility and integrity data.

Restart dumps were omitted from the published archive after the field snapshots and logs were verified. The archived results are for analysis and visualization, not direct checkpoint restart.

The temporary 6 mW restart and failed-state dump copies were removed after publication to recover 257 MB. Their paths, sizes, and SHA-256 values are recorded in `results/full_6mW_3x_safe_20261106/diagnostics/checkpoint_cleanup_6mW_20261106.csv`.
