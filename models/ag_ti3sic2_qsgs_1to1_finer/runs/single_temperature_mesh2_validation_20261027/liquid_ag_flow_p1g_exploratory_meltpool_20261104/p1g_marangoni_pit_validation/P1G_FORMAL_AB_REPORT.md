# P1G Formal A/B Heating, Marangoni Flow, and Solidification

This report records the completed formal A/B runs. `P1G_FINAL_REPORT.md` remains the earlier P1G-1 program audit and short-test record; its statement that the formal runs had not started is historical.

## Outcome

Both cases completed from the same 300 K state, reached the same peak-pool time, shut off heat at the measured Ti₃SiC₂ temperature threshold, and cooled until the liquid Ag volume was zero for two consecutive steps. Case B generated nonzero motion with Marangoni enabled, but the audited exposed surface did not form a resolved depression at any saved field frame. The maximum measured change was an upward shift of one grid cell; no permanent pit formed in this exploratory run.

The A/B comparison is not a quantitative prediction for a real Ag–Ti₃SiC₂ contact. The continuous-medium flow parameters, contact-line physics, and nanometre-scale internal flow are not sufficiently calibrated or resolved for that interpretation.

The 53 output times match between A and B. Across those matched times, the maximum absolute differences are 0.448 K in the Ti₃SiC₂ peak temperature, 5.44 K in the Ag peak temperature, and 0.00980 nm³ in total liquid-Ag volume; both cases have the same first-melt, heat-off, and solidification timestamps.

## Matched setup

- Geometry: the common Ag70 QSGS material masks, mapped from 64×64×50 to a frozen uniform 64×64×64 domain. The geometry audit records the remapping and its small dimensional changes.
- Grid spacing: 0.3018574846 nm; cubic domain edge: 19.318879 nm.
- Initial temperature: 300 K in both cases.
- Heat source: fixed 8.93189911185 mW Gaussian/depth-exponential source applied to the Ag + Ti₃SiC₂ composite; no predictive throttling.
- A: Marangoni disabled. B: Marangoni enabled with the negative-slope, reducing-atmosphere candidate.
- Heating stops at an actual Ti₃SiC₂ maximum of 2100 K or after full-power input and the ideal bottom sink differ by ≤0.01 mW for ten consecutive solver steps. The 2300 K property limit remains a hard ceiling. Both runs stopped by the temperature condition.
- Boundaries and model: x/y periodic, ideal bottom heat sink, fixed Ti₃SiC₂ skeleton, CLSVOF macro free surface. No evaporation, recoil pressure, ejection, or Ag mass-loss mechanism.
- Each case has 53 solver-history rows and 19 saved 64³ field frames from initial state through solidification.

The candidate surface-tension relation in this run is

`γ(T) = [984 − 0.298 (clamp(T, 1219, 1734) − 1234)] × 10⁻³ N/m`.

The viscosity uses the NIST Arrhenius fit, clamped to its experimental 1273–1673 K interval. Above that interval the code holds the endpoint value; the simulation therefore does not validate liquid-Ag viscosity or surface-tension behaviour at the peak temperatures.

## A/B results

| Measure | A: Marangoni off | B: Marangoni on |
|---|---:|---:|
| First liquid Ag | 0.132666250 ps | 0.132666250 ps |
| Heat off | 0.624017255 ps; scaffold 2117.734 K | 0.624017255 ps; scaffold 2117.286 K |
| Ag temperature at heat off | 2093.287 K | 2093.477 K |
| Total liquid Ag at heat off / maximum | 121.852732 nm³ | 121.854858 nm³ |
| Largest connected liquid cluster | 110.090850 nm³, 4035 cells | 110.083468 nm³, 4035 cells |
| Largest-cluster bounding dimensions | 10.263 × 10.565 × 3.320 nm | 10.263 × 10.565 × 3.320 nm |
| Maximum solver-domain speed | 0 m/s | 62.2876 m/s at 0.781874 ps |
| Fully solidified | 0.861280016 ps | 0.861280016 ps |
| Ag-volume drift, maximum absolute | 0% | 4.7105×10⁻⁵% |
| Maximum `|energy residual| / enthalpy` | 2.8720×10⁻⁸ | 2.8720×10⁻⁸ |

The B peak speed occurred during late cooling, when only 8.78 nm³ of Ag remained liquid. It is a whole-domain maximum, not a liquid-phase average. Its size should not be interpreted as a calibrated physical Ag velocity.

## Free-surface and pit audit

The earlier postprocessor calculated a volume-equivalent column height by summing VOF through z. That does not measure surface elevation for columns containing internal voids. The postprocessor has been corrected to use the highest z cell with `VOF > 0.5` in each x-y column, matching the solver’s top-envelope height diagnostic. A synthetic test verifies the highest-cell and surface-temperature selection, and both HDF5 files and pit tables were regenerated.

Pit/rim classification uses a displacement threshold of 0.5 grid spacing (0.150929 nm). In Case B, every saved frame has zero pit depth, zero pit area, and zero pit volume above that threshold. The largest upper-envelope uplift is 0.301857 nm (one grid cell); the source-region mean elevation change is +0.173340 nm at heat-off and +0.202756 nm after solidification (positive means upward). Case A remains stationary. Thus the computed result is slight surface uplift at this resolution, not a central depression or a permanent pit. Sub-grid surface motion cannot be distinguished reliably on this mesh.

The upper-envelope audit applies only to the exposed outer sample surface, assumed to be single-valued in z. It does not resolve Ag/gas interfaces inside the Ti₃SiC₂ pore network, overhangs, or individual nanoscale channels.

## Energy and mass checks

At final solidification, Case A records 5.732795×10⁻¹⁵ J cumulative input, 1.293280×10⁻¹⁵ J ideal-sink removal, and 4.439515×10⁻¹⁵ J stored enthalpy. Case B records 5.732818×10⁻¹⁵ J input, 1.293279×10⁻¹⁵ J sink removal, and 4.439539×10⁻¹⁵ J stored enthalpy. The maximum absolute residuals in the histories are 8.41×10⁻²⁴ J (A) and 3.55×10⁻²⁴ J (B). Ag volume is conserved closely; B’s maximum absolute drift is 4.71×10⁻⁵%.

The bottom boundary is an ideal heat sink. The cooling phase represents internal redistribution and this idealized sink, not a calibrated real electrode or ambient heat-loss model.

## Numerical and physical limitations

- The source geometry maps the whole composite to one VOF sample/gas interface. Internal pore gas/liquid surfaces are absent; this run cannot answer how Ag moves through real pores.
- The frozen 64³ mesh has 0.302 nm cells. A previously measured 1.55 nm neck spans about five cells, insufficient for a quantitative pore-flow claim. No grid-convergence or gas-cap-height sensitivity study was done for this coupled case.
- Surface-tension and viscosity data are candidates with clamped temperature ranges, not validated at the run’s peak Ag temperature. The contact angle and dynamic contact-line friction are also not calibrated here.
- The logged pressure residual is not nondimensionalized or tied to a validated acceptance threshold. The successful process exit and small energy/mass errors do not by themselves validate the velocity magnitude.
- Basilisk `qcc` dimensional checking was disabled because its dimension parser aborted; units were manually audited and HDF5 attributes were checked.
- The run is an exploratory isolated prototype. It did not alter the production enthalpy solver, and it does not establish grid convergence or a real-material pit prediction.

## Reproducibility and outputs

- Source: `p1g_free_surface_coupled_audited.c`; restart support: `two-phase-clsvof-p1g.h`.
- Basilisk snapshot: `comphy-lab/basilisk-C v2026-07-20`, commit `9bb4f30aaf195feaa7a3a613468a1c888f5a4b9e`.
- Both formal case `code_hashes.sha256` files match for the solver source, restart header, executable, and input manifest.
- `results/formal_case_A_.../` and `results/formal_case_B_.../` contain run configuration, logs, 19 compressed field snapshots, history and audit CSVs, checkpoint sidecars, heat-off/solidification dumps, full HDF5 fields, and maximum-pool 3-D images. The full HDF5 and 70 MB restart dumps remain in the cloud workspace. `thermal_flow_fields_selected.h5` contains the 19-frame thermal, phase, geometry, pressure, velocity, candidate surface-tension, and candidate Marangoni-acceleration fields for sharing; it omits derived vorticity/interface normals and the full enthalpy-density field. Its sizes are 15.7 MB (A) and 72.8 MB (B), below GitHub’s 100 MB individual-file limit.
- `results/P1G_A_B_COMPARISON.csv` and `results/P1G_A_B_COMPARISON.png` provide the matched time-series comparison.
- `results/P1G_PIT_GEOMETRY.csv` in each case uses the corrected upper-envelope method and a half-cell pit/rim threshold.
- `results/formal_case_A_.../maximum_pool_3d.png` and `results/formal_case_B_.../maximum_pool_3d.png` show the 3-D phases at their maximum total liquid-Ag frame.

To regenerate derived fields and the comparison plot in the run environment:

```bash
PYTHONPATH=/tmp/p1g-python-deps python3 p1g_postprocess.py results/formal_case_A_marangoni_off_2100K_cut_cool_20261105 --output results/formal_case_A_marangoni_off_2100K_cut_cool_20261105/thermal_flow_fields.h5
PYTHONPATH=/tmp/p1g-python-deps python3 p1g_postprocess.py results/formal_case_B_marangoni_on_2100K_cut_cool_20261105 --output results/formal_case_B_marangoni_on_2100K_cut_cool_20261105/thermal_flow_fields.h5
PYTHONPATH=/tmp/p1g-python-deps MPLCONFIGDIR=/tmp/p1g-mpl python3 compare_p1g_ab.py .
PYTHONPATH=/tmp/p1g-python-deps python3 test_p1g_postprocess_surface.py
```
