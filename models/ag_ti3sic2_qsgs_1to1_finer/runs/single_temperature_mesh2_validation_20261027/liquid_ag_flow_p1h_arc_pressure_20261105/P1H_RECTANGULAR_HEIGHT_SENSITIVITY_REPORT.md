# P1H rectangular-domain height and top-boundary sensitivity

Date: 2026-11-05

Branch baseline: `codex/ag70-p1h-arc-pressure` at `f62aa755`
Purpose: determine whether the short Ag70 free-surface response changes when only the gas-cap height or top fluid boundary changes.

## Setup

The new isolated solver, `p1h_free_surface_rectangular.c`, uses Basilisk's fixed `multigrid3D` grid. It reads the existing 64³ Ag70 material/source maps without modifying them, keeps the x/y periodic box at 19.318879 nm and the cell spacing at 0.301857485 nm, and fills added z cells with gas and zero source. The unchanged sample occupies 53 z cells (Ag 151,809 cells; Ti₃SiC₂ 65,279 cells; sample volume 5,914.596168 nm³).

The six matched runs used 64×64×64, 64×64×76 and 64×64×86 cells, each with closed and ambient-pressure open top boundaries. Their gas caps are 11, 23 and 33 cells (3.3204, 6.9427 and 9.9613 nm). The middle cap is one cell taller than the originally proposed 22-cell cap because 64×64×75 has no common dyadic coarsening level in this multigrid layout. The exact 75-cell attempt was preserved and failed: after the first step, the temperature solve reached residual `1.20546e181` at 200 iterations and the process terminated with a floating-point exception. The adjacent 76-cell grid has two coarsening levels and completed. The exact 86-cell grid has one level and completed, but required substantially more runtime than 76 cells.

The thermal source and boundary data match the existing 0.16 ps software smoke: 8.931899 mW source power with a 0.02 ps ramp, Marangoni disabled, and a uniform 1.8 MPa surface-pressure code-path test with the same ramp. The 1.8 MPa value is a numerical smoke input only; it is not calibrated as a real nanoscale arc pressure. Each case stopped at 0.16 ps. No cooling stage or physical crater prediction was attempted.

## Results

The table reports the last history row at 0.14778175 ps; the complete field snapshots and restart checkpoints were written at exactly 0.16 ps. The liquid-volume number is the integrated Ag liquid-fraction volume.

| Grid | Gas cap | Top | Peak T (K) | Peak Ag (K) | Peak skeleton (K) | Liquid Ag (nm³) | Max effective pressure (Pa) | Integrated z force (N) | Max speed (m/s) | Depression (nm) |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 64×64×64 | 11 cells / 3.320 nm | closed | 1287.09664 | 1244.79318 | 1287.09664 | 0.190788090339 | 718992.111665 | −1.27422663929e−13 | 0.003078114 | 0 |
| 64×64×64 | 11 cells / 3.320 nm | open | 1287.09664 | 1244.79318 | 1287.09664 | 0.190788090337 | 718992.111654 | −1.27422663927e−13 | 0.003078114 | 0 |
| 64×64×76 | 23 cells / 6.943 nm | closed | 1287.09373 | 1244.79113 | 1287.09373 | 0.190770978162 | 718878.709069 | −1.27402239990e−13 | 0.003084051 | 0 |
| 64×64×76 | 23 cells / 6.943 nm | open | 1287.09373 | 1244.79113 | 1287.09373 | 0.190770978146 | 718878.708954 | −1.27402239968e−13 | 0.003084051 | 0 |
| 64×64×86 | 33 cells / 9.961 nm | closed | 1287.09373 | 1244.79113 | 1287.09373 | 0.190770940764 | 718878.511959 | −1.27402201761e−13 | 0.003119417 | 0 |
| 64×64×86 | 33 cells / 9.961 nm | open | 1287.09373 | 1244.79113 | 1287.09373 | 0.190770945456 | 718878.540346 | −1.27402207195e−13 | 0.003119417 | 0 |

Across 11–33 gas-cap cells, peak temperature changes by 0.00291 K and liquid volume by about 0.0090%. The maximum effective pressure and integrated pressure force change by about 0.016%. At each height, open and closed top results agree to the reported precision: liquid-volume differences are below 1.6×10⁻⁸ nm³ and no surface depression is resolved. Maximum Ag-volume drift over all time histories is 1.78×10⁻¹⁴%; the maximum absolute energy residual is 2.76×10⁻²² J. All six final temperature/VOF/material/pressure snapshots passed binary-size, dimension, finite-value, and fraction-range checks.

The fixed 64³ multigrid run also reproduces the existing 64³ octree smoke closely: peak temperatures agree to the printed precision, liquid volume differs by 4.4×10⁻¹⁰ nm³, and the effective pressure differs by 0.0026 Pa. This is a useful baseline check for the new rectangular solver.

## Decision and limits

The z-extended domain works without changing the sample geometry, x/y period or voxel spacing. For this early heating window, neither the 11-cell cap nor the top boundary measurably affects the sample response. The 22-cell target is represented by 23 cells to retain a stable multigrid hierarchy; the 33-cell target is exact.

The runs show only the first small liquid cluster near the Ag melting range and no resolved surface depression. They do not answer whether Marangoni flow forms a permanent crater: Marangoni was off, the surface-pressure value was only a code-path smoke input, and no cooling/solidification was run. The next physical exploration still needs a separately justified pressure input (or a no-pressure Marangoni-only control), followed by a longer heating and cooling run on the validated rectangular geometry.

## Files and reproduction

- `run_p1h_rectangular_height_sensitivity.sh` builds and runs the six cases with the pinned Basilisk source.
- `results/rectangular_height_sensitivity_20261105/rectangular_height_sensitivity_summary.csv` contains the time-series summary.
- `results/rectangular_height_sensitivity_20261105/snapshot_validation.csv` records final-snapshot shape and integrity checks plus SHA-256 hashes.
- Each `nz{64,76,86}_{closed,open}/` folder contains its run config, histories, logs, compressed 0.16 ps field/pressure frames, and compressed 0.16 ps restart dump with completion/state sidecars. The tracked result set includes the final field/pressure frames and end-window checkpoints; it omits initial/intermediate frames and checkpoints plus the duplicate `p1h_final.dump.gz`. `latest_checkpoint.txt` points to the recoverable restart dump. Before resuming one locally, decompress the pointed-to file with `gzip -dk`.
- `results/rectangular_geometry_smoke_nz75_20261105/run.stderr.log` preserves the exact 75-cell failure.

Reproduction command (from this directory):

```bash
BASILISK_SRC=/path/to/pinned/basilisk/src \
OMP_NUM_THREADS=5 \
bash ./run_p1h_rectangular_height_sensitivity.sh
```

The run is intentionally a short 0.16 ps preflight, not a production heating/cooling result.
