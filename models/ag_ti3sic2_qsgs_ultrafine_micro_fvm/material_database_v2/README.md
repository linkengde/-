# Ag/Ti₃SiC₂ material database v2

This folder contains an independently versioned, literature-audited property set and its FVM checks. It does not overwrite the earlier material tables, geometry, MD inputs, or potential files. Supplied papers are treated as scientific sources; no paper PDFs are copied here.

## Reproduce the property checks

Install the optional Numba dependency with the project requirements, then expose the target directory when running Python:

```bash
cd /workspace/-/models/ag_ti3sic2_qsgs_ultrafine_micro_fvm/solver_optimization
python -m pip install --target /tmp/solver-opt-deps --upgrade -r requirements-numba.txt
cd ../..
PYTHONPATH=/tmp/solver-opt-deps python -m unittest material_database_v2.test_material_database -v
```

## Main run scripts

- `build_database.py` rebuilds the CSV tables from the recorded source equations and checkpoints.
- `compare_50_old_vs_pankratz.py` runs the matched 50³ old/new material heat cycle.
- `run_material_validation_suite.py` checks Stefan fronts, energy balance, and interface resistance.
- `run_width_sensitivity_50.py` compares 0.5, 1, and 2 K numerical mushy widths.
- `run_100_time_refinement.py` runs the 100³ Pankratz v2 heat/cool cycle with chosen time-step factors. It writes restart checkpoints after accepted steps. Start with factors `1 0.5`; run `0.25` only if the half-step metrics differ by more than 3%.
- `run_n100_freezing_event_refinement.py` restarts the 100³ cooling from its saved heat-end field and caps the step while any Ag remains liquid; use a 5/10 ps pair to resolve the complete-freeze event without rerunning heating.
- `run_150_short_melt_probe.py` stops at the first melt or at a strict step/wall cap. It does not start a full 150³ cycle.
- `run_scaled_phase_case.py --label micro200um --edge-scale 20` runs the 200 μm scale-similarity case: source density scales as `L⁻²` and both physical durations/time steps scale as `L²`. Use `--source-scale` only for the separately labelled energy-budget probe.
- `run_scaled_phase_case.py --label exploratory_target1400K --edge-scale 1 --source-scale 1.176` runs the independent 50³ numerical energy-budget pulse. It is not an arc load.
- `analyze_source_region.py` reports source-region and full-RVE Ag liquid fractions separately, plus 26-neighbor connected regions above explicit Ag/liquid thresholds.
- `compare_scale_similarity.py` writes the field-difference and event-time scale-ratio record for the 10/200 μm pair.
- `compare_freezing_event_refinement.py` summarizes the 10/5 ps cooling-event pair alongside the uniform-step cycles.
- `render_material_figures.py` regenerates literature and computed-field figures.

For the 100³ and 150³ runs, use the already installed Stage 4 Numba/EW-safe/FFT solver, four worker threads, and numerical validation pulse documented in the result JSON files. The pulse is not a measured arc boundary. Liquid flow, evaporation, recoil pressure, and electrical heating are outside this phase-change validation.

## Evidence and limitations

`MATERIAL_PROPERTY_AUDIT.md` records source ranges, units, conversions, and discrepancies. `PROPERTY_COMPARISON.md` records same-load numerical results and acceptance decisions. `SOURCE_MANIFEST.json` identifies the supplied literature without redistributing PDFs. Raw accepted-step CSVs, NPZ fields, and restart checkpoints are under `results/`; figures are under `figures/`.

The Pankratz set is the active thermodynamic baseline. Feber and Cagran remain separate alternatives. The liquid-Ag conductivity is Ho et al.'s provisional estimate with its stated uncertainty. Ag/Ti₃SiC₂ interface resistance has no direct measurement among the supplied sources; nonzero values are numerical sensitivity cases only.

The 50³ old/new same-load cycle and the independent stronger numerical pulse both fully cool and freeze. The stronger pulse yields one connected region of 51 cells above Ag fraction 0.5 and liquid fraction 0.5; this remains a solver validation, not an experimental melt pool. All three 100³ uniform time-step cycles complete, but complete-freeze time is not converged within 3% across adjacent uniform schedules. A 5/10 ps event-focused cooling pair agrees within 0.98% only when started from the same nominal-step heated field. The 150³ run stopped at its first mushy cell and is not evidence for a full cycle. See the decision register in `PROPERTY_COMPARISON.md` before using these results.
