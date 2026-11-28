# Property-set comparison

| Property | Active v2 main set | Separate literature comparison | Use in solver |
|---|---|---|---|
| Ag melting point | Pankratz, 1235.08 K | Abdullaev, 1234.93 K | Pankratz value |
| Ag solid Cp | Pankratz polynomial, 298.15–1235.08 K | Cagran high-T fit, 1000–1234 K | Pankratz polynomial |
| Ag liquid Cp | Pankratz, 33.472 J/(mol K) = 310.30 J/(kg K) | Feber, 314.61 J/(kg K); Cagran, 259.8 J/(kg K) | Pankratz; Feber/Cagran remain independent sensitivities |
| Ag heat of fusion | Pankratz, 11296.8 J/mol = 104.728 kJ/kg | Feber: 2590 or 2710 cal/mol; Cagran: about 101.6 kJ/kg text/Table 2 | Pankratz; no averaging |
| Ag solid k at melt | Ho recommended, 358 W/(m K) | Old candidate incorrectly used 175 W/(m K) | 358 W/(m K) |
| Ag liquid k at melt | Ho provisional estimate, 175 W/(m K), about 20% uncertainty | Kondo cites about 177 W/(m K) near 1273 K via Wiedemann–Franz | 175 W/(m K); sensitivity advised |
| Ag density | Fixed 10499.4 kg/m³ at 300 K | Abdullaev temperature-dependent solid/liquid fits | Fixed in this solver; variable data stored for a later conservative-flow formulation |
| Ag electrical resistivity | Not coupled | Matula recommended table; Cagran high-T fits | Database only |
| Ti₃SiC₂ Cp and k | Barsoum fits, with source range limits | Existing user CSV extends fits to 2300 K | Solid skeleton; solver stops above 1473.15 K |
| Ag/Ti₃SiC₂ interface R | 0 m² K/W baseline | 1e−10, 1e−9, 1e−8 m² K/W numerical sensitivity cases | Not an experimentally measured property |

## What changed numerically in the first matched 50³ cycle

The QGS seed 1105 field, 10 μm edge, 70/30 phase fraction, 4e18 W/m³ localized numerical pulse, 1 ns heating, 2 μs cooling, 0.25 nominal time-step factor, initial/boundary 300 K, zero interface resistance, nonlinear/linear tolerances 1e−8, and Numba/EW-safe/FFT settings were held constant. These are validation loads, not an arc history.

| Metric | Old candidate | Pankratz v2 | Interpretation |
|---|---:|---:|---|
| Peak temperature | 1242.611 K | 1235.357 K | Correcting solid Ag k reduces the hot spot above the melting point |
| First nonzero Ag liquid fraction | 0.725 ns | 0.800 ns | Discrete first-crossing time at the common 25 ps heating step |
| Maximum cell Ag liquid fraction at heat end | 1.000 | 0.7766 | Local cell state |
| Ag bulk liquid fraction at heat end | 3.7846e−4 | 1.5517e−4 | Volume-weighted over Ag only |
| Ag liquid volume at heat end | 2.6492e−19 m³ | 1.0862e−19 m³ | Numerical QGS RVE volume |
| Ag / Ti₃SiC₂ mean temperature at heat end | 303.907 / 300.819 K | 304.032 / 300.892 K | Includes the complete phase masks |
| Latent energy at heat end | 2.9105e−10 J | 1.1930e−10 J | Liquid-volume estimate using fixed reference Ag density |
| Complete-freeze crossing after cooling starts | 25 ns | 25 ns | First cooling output with Ag bulk liquid fraction ≤ 1e−8 |
| Final Ag bulk liquid fraction | 0 | 0 | Both runs finish fully frozen |
| Final normalized energy residual | 9.98e−10 | 2.35e−9 | Both are far below the 0.1% acceptance limit |

Both cases complete 120 accepted steps. The matched load produces only a very small partial liquid volume under v2; it does not establish a meaningful melt pool. No input was adjusted to restore the old melt fraction. A separate, explicitly numerical energy-budget study is required before making a stronger phase-change validation claim.

Full per-step CSVs, restart checkpoints, final fields, heated-endpoint fields, and run JSON are in `results/same_load_50/`. This comparison isolates the material change only for the specified 50³ numerical cycle; it is not a real contact or arc-performance prediction.

## Stage 5 validation results

### Property model and small verification suite

The v2 CSV/JSON tables and Python model preserve the Pankratz main set and keep Feber and Cagran as separately sourced alternatives. Cagran's liquid heat capacity of 259.8 J kg⁻¹ K⁻¹ is retained with the paper's stated 1234–2000 K fit range. Its solid enthalpy/table mismatch is recorded rather than silently reconciled. Ho's solid Ag conductivity is now about 358 W m⁻¹ K⁻¹ at melting; the prior candidate's 175 W m⁻¹ K⁻¹ was incorrectly applied to solid Ag. Liquid conductivity remains a provisional Ho estimate with about 20% uncertainty. No measured Ag/Ti₃SiC₂ interface resistance was found.

The 9 Python verification tests pass: enthalpy/temperature round trips, Cp–enthalpy derivative, latent heat counted once, range-limit errors, Numba/Python field agreement, one-cell energy and 300 K hold, and checkpoint restart identity. The validation suite reports a fine-grid Stefan melt-front error of 0.815% and freeze-front error of 5.414%; the former is within a 2% target, the latter is not. The qR interface-jump numerical check agrees to relative error below 4e−12 for the nonzero R probes. All liquid fractions stay in [0, 1], and the tested normalized energy residuals are orders of magnitude below the 0.1% criterion. Thus conservation and implementation checks pass, while the Stefan freezing-front accuracy remains a validation limitation.

The 0.5, 1, and 2 K numerical mushy widths were retained as distinct runs. Heating-end bulk Ag liquid fractions span 1.5492e−4 to 1.5567e−4, and peak-temperature changes are below 0.02%; no width change was used to accelerate the solver.

### Same-load 50³ comparison and independent melt probe

The matched old/new 50³ cycle completes heating and cooling under identical QGS seed 1105, source, time schedule, boundary and solver tolerances. Correcting solid Ag k raises its melt-point value from 175.14 to 358.04 W m⁻¹ K⁻¹. The v2 hot spot reaches 1235.357 K, with 0.0001552 full-RVE Ag liquid fraction and 0.05281 source-region Ag liquid fraction; the old candidate gives 0.0003785 and 0.1288, respectively. Both freeze fully and have maximum normalized energy residual below 3e−9. The unchanged-load result is intentionally reported as a small partial melt, not a developed melt pool.

An independently labelled 50³ energy-budget pulse at 4.704e18 W m⁻³, integrated input 9.676e−9 J, produces a maximum temperature of 1393.405 K. At heat end, Ag liquid fraction is 0.0005929 over the full RVE and 0.2018 in the source cells, with 4.150e−19 m³ Ag liquid volume. At thresholds Ag fraction ≥0.5 and liquid fraction ≥0.5, 51 connected cells form one 26-neighbor component (4.08e−19 m³ geometric cell volume; Ag-weighted liquid volume 3.641e−19 m³). The run completes full cooling, returns to zero Ag liquid fraction, and has maximum normalized energy residual 8.01e−10. This establishes a useful numerical phase-change probe, not an experimental load or arc simulation.

### 100³ cycles, time refinement, and scale similarity

The 100³, 10 μm cases at nominal, half and quarter time-step factors all complete a 1 ns numerical heating pulse followed by 2 μs cooling. Peak temperature is 1235.383–1235.422 K, heat-end bulk Ag liquid fraction is 1.2145e−4–1.3339e−4, and each run fully freezes. Peak RSS is 520–564 MiB and wall time is 12.1–13.8 min. Maximum normalized cumulative energy residual is at most 2.86e−9.

Uniform fixed-step refinement does not satisfy the requested 3% criterion for all integrated phase metrics: from factor 0.5 to 0.25, heat-end bulk liquid changes 2.87%, but cooling complete-freeze time changes 8.02%; from factor 1 to 0.5 those changes are 6.76% and 19.95%. The 50³ and 100³ source masks also represent the localized pulse differently: integrated source energy is 8.228e−9 J and 8.020e−9 J, a 2.52% geometry/discretization difference. These effects must not be conflated as temporal error or claimed as spatial convergence.

For a local event check, 10 ps and 5 ps maximum-step cooling reruns start from the same saved nominal-step heated 100³ field and carry forward the initial enthalpy and source-energy bookkeeping. Complete-freeze time differs by 0.98% between those two event-capped cooling runs. This supports local cooling-event resolution only; it does not remove the mismatch with uniform-step full-cycle runs, nor validate the event time for a quarter-step heated field.

At 150³, a bounded first-onset probe reaches the numerical mushy interval after 0.776 ns, 17 accepted steps and 265 s, with about 1.125 GiB peak RSS. It has only one incipient mixed cell (maximum liquid fraction 0.00110, full-RVE Ag liquid fraction 4.65e−10, liquid volume 3.26e−25 m³), not a meaningful connected liquid region. A full 150³ heat/cool cycle was not started.

The 10 μm and 200 μm 50³ QGS models provide a strict scale-similarity check when source density and time scale as L⁻² and L². Heated-end temperature fields differ by at most 4.55e−13 K; phase masks and liquid fields are identical. First-melt and complete-freeze times are exactly 400× in the recorded outputs, matching (20)². This validates implementation of the similarity transformation only; it does not compare real arc loading or contact performance.

### Decision register

| Item | Decision | Evidence / limitation |
|---|---|---|
| Literature-backed v2 database | GO with stated ranges | Main Pankratz set; Feber/Cagran remain separate alternatives; source discrepancies and extrapolations are explicit. |
| Old solid Ag conductivity correction | GO | 175 → about 358 W m⁻¹ K⁻¹ at melting, with separate liquid estimate. |
| Material model, conservation, checkpoints | GO for tested cases | 9 tests pass; energy residuals <3e−9; interface qR check passes. |
| 50³ same-load melt/cool | Conditional GO | Full cycle and freeze pass, but v2 produces only a small partial liquid region. |
| 50³ stronger numerical melt probe | GO as a code/physics-path test | One connected thresholded region of 51 cells; not a real arc or validated experimental melt pool. |
| 100³ complete heat/cool execution | GO for the tested numerical pulse | All 3 time-step runs complete with low residuals and measured resource use. |
| 100³ time convergence at 3% | NO-GO | Uniform schedules miss the target for cooling freeze time and, in some pairs, integrated liquid. Local 5/10 ps event pair agrees within 0.98% from the same nominal heated field only. |
| 100³ spatial/grid convergence | NOT ESTABLISHED | Source-mask energy changes 2.52% between 50³ and 100³; local neck and melt-boundary convergence were not assessed. |
| 150³ first mushy-onset probe | Conditional GO | Bounded probe completed; detected only a tiny incipient mixed-cell liquid amount. |
| 150³ full cycle | NO-GO for this stage | Onset probe already takes 265 s/1.125 GiB; no demonstrated full-cycle performance budget. |
| Real arc/contact prediction, liquid flow, Marangoni, evaporation | NO-GO / out of scope | No measured local heat history, direct interface R, or validated liquid-flow model in this stage. |

The computed outputs and plotting utility are under `results/`, `figures/`, and the `run_*.py` scripts in this folder. `LOCAL_SOURCE_REGION_METRICS.json` separates local source-cell values from whole-RVE liquid fractions; `SCALE_SIMILARITY_COMPARISON.json` records the dimensionless-scale check. All applied source terms are numerical validation pulses.
