# V9–V16 cumulative model review and V17 recommendations

**Scope.** This is a retrospective comparison of the archived V9–V16 checkpoints. All eight models are scored, or have hash-checked archived scores reused, on the same six already-scored registry references: the three V15 references and the three V16 references. No new labels were read, and no DFT, training, MD, TTM, recalibration, or model selection was performed.

## Findings for A

1. **There is real improvement on the matched six-reference set, but none of V9–V16 meets the per-interface force-vector gate.** Pooled force-vector RMSE falls from 0.2905 eV/Å (V9) to 0.1464 (V16), a 49.6% reduction. V16 is also 21.0% below V15 on this same set. V16 per-interface pooled RMSE is Ag–C 0.1048, Ag–Si 0.2096, and Ag–Ti 0.0905 eV/Å; all exceed 0.05.
2. **Energy error improves sharply after reference calibration, but this does not explain force error.** V10–V16 pooled energy MAE on the six references falls from 211.0 to 2.20 meV/atom. V9 energy is not interpretable as a calibrated comparison. The V16 own-holdout energy gates pass, while all three force gates fail.
3. **The residual is localized and directional.** On the V16 Ag–C holdout, Ag#688433 has a 0.445 eV/Å vector error, mostly transverse (0.431 eV/Å); nearby framework Ti#604805 is 0.291 eV/Å. On the V16 Ag–Ti holdout, marked Ti#479027 is 0.317 eV/Å and contributes −0.282 eV/Å to the longitudinal separation error. On the V15 Ag–Si reference, the V16 checkpoint still has a 1.135 eV/Å transverse error on marked Si#366934 and 1.344 eV/Å on nearby Ag#685172, although its scalar separation error is only 0.0079 eV/Å. Therefore a passing separation component can conceal a large vector residual.
4. **The strongest supported V17 choice is to keep the V16 recipe fixed while using the eight registered signed diagnostic labels.** They can test the targeted local response around four neighborhoods. They cannot establish independent generalization, resolve inherited label provenance, or justify MD/TTM. Do not increase force weight or epochs without a fixed-data, one-factor comparison.

## Provenance and comparison method

The driver is `analyze_v9_v16_common.py`. The machine-readable results are under `analysis_04/`; the six reference structure and input SHA256 values, model SHA256 values, scoring runtime, metrics, and source evaluation hashes are recorded in `analysis_04/common_geometry_scores.json`. The flat `analysis_04/common_geometry_scores.csv` has one row for each of 8 models × 6 references (48 rows). `analysis_04/version_ledger.csv` and `.json` preserve the version-by-version source paths, hashes, dataset counts, known training settings, and evaluation-role limitations.

Scoring used CPU float64, one Torch thread, PyTorch 2.5.1+cpu, mace-torch 0.3.16, and ASE 3.29.0. Force-vector RMSE is `sqrt(mean_i(sum_xyz((Fpred-Fref)^2)))`. Marked-pair separation is the relative force projected along the direct marked Ag–X vector used by the published scorer. “Local” is the existing MIC distance mask within 4.5 Å of either marked atom.

V15 and V16 scores were reused only where the prior matched run had the same model hash, reference hash, runtime convention, and metric; their published scores were reproduced. Three prior V14/V15 comparisons on V15 references were likewise reused after hash checks; V14 was inferred on the V16 references. V9–V13 values on these six references are new retrospective scores, not reproductions of those versions’ original evaluations. No holdout role has been relabeled: the V15 and V16 structures retain their original role-at-first-scoring in the CSV and JSON.

The shared set is useful for a controlled comparison of frozen checkpoints, but it is **not** a fresh independent test set. The V15 geometries are small perturbations of prior training-parent motifs; V16 references are local-registry structures; the old per-version training membership is not completely reconstructable; and the V14 Ag–Ti holdout exactly duplicates a V14 training geometry. These data do not test other interface morphologies, thermal disorder, liquid Ag, pressure, or extended-cell transfer.

## V9–V16 data and recipe timeline

Counts below are `(train, valid, test)` and label counts refer to split records. Unknown settings remain unknown; version names are not evidence of an isolated experimental factor.

| Version | Dataset and evidenced change | Training/selection evidence | Evaluation caveat |
|---|---|---|---|
| V9 | Initial four-element, force-focused model. Dataset, splits, label counts, and complete launch command are absent from the handoff. | Checkpoint/hash available; recipe and selected epoch unknown. | Original membership is not fully reconstructed; energy calibration was not established. |
| V10 | Packaged interface set is 23/2/3, with 15 train energy labels and 23 force labels. Added four PBE elemental anchors, Ag–Ti 2.55 Å energy/force, and an energy label for the 2.3871 Å frame; retained the 2.70 Å same-motif holdout. Energy-composition rank is 4/4. | Exact launch settings and per-version log are unavailable. | Elemental reference meshes differ from interface meshes; old model membership is incomplete. |
| V11 | Checkpoint is available; its data directory documents the V10 23/2/3 dataset (15/23 train energy/force labels). | “Energy-focused” is documented, but no isolated V10→V11 recipe change or separate launch log is available. | Cannot attribute V10/V11 score movement to a specific parameter. |
| V12 | 23/2/3; train energy/force labels 19/23. Replaced four Ag–Si/Ag–C force-only LCAO frames with same-geometry PW-PBE energy/force labels. | MACE-MP-0b3 medium, seed 45, CPU, batch 1, max 80 epochs, lr 1e-4, wd 5e-7, energy weight 100, force weight 1000, estimated E0s; selected epoch 76 as recorded in history. | Three frozen tests and a small correlated external set; historical training metric included a force-only energy-mask issue. |
| V13 | 26/2/3; train energy/force labels 22/26. Added three verified PW-PBE distance labels: AgC 2.30 Å, AgSi 2.60 Å, AgTi 2.60 Å. Existing valid/test bytes retained; composition rank 4/4. | Same documented V12 recipe; 80 epochs, selected epoch 72. | Some midpoint checks share parent motifs with added training distances; local-registry screens still failed force gates. |
| V14 | 31/2/3; train energy/force labels 27/31. Added five PW-PBE labels: two lateral-registry structures, an AgTi probe, and AgC/AgSi strain structures. | Same documented recipe; selected epoch 71 from the training record. | The AgTi “holdout” exactly duplicates an AgTi training geometry. |
| V15 | 29/2/3; all train rows have energy and force labels. Excluded four provenance-unknown force-only LCAO rows and one unresolved AgTi force-only row; added three verified residual-shell PW-PBE labels. | Same recipe; 80 epochs; selected epoch 75 (zero-based 74); training exit 0 and selection independent of blind labels. | Its three registry references are near parent motifs, not independent morphology tests. Provenance for inherited labels remains partly unresolved. |
| V16 | 35/2/3; all train rows have energy and force labels. Added three wide-span registry and three wide-span distance PW-PBE labels; metadata counts were repaired without changing split bytes. | Same recipe; 80 epochs; selected epoch 80 (zero-based 79). Full log/checkpoint are verified; launcher exit receipt is unavailable and is not inferred. | Three local-registry references are narrow; six-reference retrospective scores do not constitute a fresh test. |

V12–V16 retained the same documented core architecture and training controls. The dominant systematic differences are data composition/coverage and checkpoint selection. Those factors changed together, so the trend is not a causal experiment.

## Matched six-reference results

### Pooled force-vector RMSE

Values are eV/Å, pooled over all 296 atoms in the six references. The gate is 0.05 eV/Å per interface; this pooled number is descriptive and is not itself a pass gate.

| Model | Pooled | Ag–C | Ag–Si | Ag–Ti |
|---|---:|---:|---:|---:|
| V9 | 0.2905 | 0.2081 | 0.4234 | 0.1609 |
| V10 | 0.2940 | 0.2037 | 0.4342 | 0.1557 |
| V11 | 0.2836 | 0.2066 | 0.4104 | 0.1608 |
| V12 | 0.2255 | 0.2061 | 0.2958 | 0.1445 |
| V13 | 0.2251 | 0.1760 | 0.3116 | 0.1469 |
| V14 | 0.1957 | 0.1711 | 0.2588 | 0.1315 |
| V15 | 0.1852 | 0.1610 | 0.2408 | 0.1335 |
| V16 | 0.1464 | 0.1048 | 0.2096 | 0.0905 |

V12 improves much of the Ag–Si force aggregate while energy calibration improves strongly. V13 reduces Ag–C force error but Ag–Si rises on this shared set. V14 improves all three force aggregates relative to V13. V15 improves Ag–C/Ag–Si modestly while Ag–Ti is essentially flat. V16 improves all three relative to V15, but every interface remains above 0.05.

### Energy and separation trends

Energy values are MAE in meV/atom, and separation values are maximum absolute marked-pair error in eV/Å across each interface’s two structures. Each sequence is ordered V9 → V16. V9 energies are shown as archived numerical output but are **not calibrated/comparable**; interpret the calibrated energy trend from V10 onward.

| Interface | Energy MAE sequence | Maximum separation-error sequence |
|---|---|---|
| Ag–C | 492.81, 298.17, 201.01, 12.05, 1.22, 1.05, 1.89, 1.11 | 0.5091, 0.5145, 0.4948, 0.3671, 0.2772, 0.3150, 0.3874, 0.1027 |
| Ag–Si | 122.06, 206.70, 100.98, 17.13, 5.21, 7.51, 4.71, 4.45 | 0.2453, 0.2264, 0.2125, 0.1668, 0.0755, 0.0446, 0.0363, 0.0793 |
| Ag–Ti | 76.82, 128.17, 13.19, 3.71, 2.99, 0.38, 1.84, 1.03 | 0.6329, 0.6306, 0.6886, 0.5061, 0.5744, 0.4632, 0.6236, 0.3227 |

Energy calibration improves sharply by V12 and stays below 10 meV/atom on this retrospective set from V13 onward. Separation does not improve monotonically: Ag–C remains marginal/poor in one V16 reference, Ag–Si is low on these structures, and Ag–Ti improves but remains the largest V16 separation error. A good scalar separation score is not evidence that all atom vectors are accurate.

### Current V16 holdout gates

This table is the V16 model’s **own three holdouts**, separate from the maximum-across-two-structures summary above. Gates are energy ≤10 meV/atom, force-vector RMSE ≤0.05 eV/Å, and separation error ≤0.10 eV/Å.

| Interface | Energy MAE (meV/atom) | Force-vector RMSE (eV/Å) | Separation error (eV/Å) | Gate result |
|---|---:|---:|---:|---|
| Ag–C | 0.698 | 0.1295 — fail | 0.0417 — pass | Force fails |
| Ag–Si | 3.021 | 0.0961 — fail | 0.0793 — pass | Force fails |
| Ag–Ti | 1.612 | 0.1053 — fail | 0.3227 — fail | Force and separation fail |

Thus two of the three V16 separation checks pass, while all three vector-force checks and Ag–Ti separation fail. The six-reference matrix is not a substitute for this own-holdout statement.

## Atom-level evidence and residual modes

- **Ag–C:** In V16’s own holdout, marked Ag#688433 has vector error 0.4451 eV/Å with 0.4305 eV/Å transverse to the marked Ag–C axis. Local RMSE is 0.1507 versus 0.0830 outside 4.5 Å. Framework Ti#604805 contributes a 0.2915 eV/Å vector error. This points to a registry/lateral response plus its adjacent shell; the low separation error (0.0417) alone hides it.
- **Ag–Si:** On the already-scored V15 Ag–Si reference, V16 has full force RMSE 0.2804 and separation error only 0.0079 eV/Å. Marked Si#366934 has 1.1349 eV/Å transverse error; nearby Ag#685172 has 1.3440 eV/Å total error. On V16’s own Ag–Si holdout, force RMSE is 0.0961 despite separation error 0.0793 passing. The transverse/interface response is the clearest persistent V16 weakness; the pair projection does not cover all of it.
- **Ag–Ti:** On V16’s own holdout, marked Ti#479027 has 0.3172 eV/Å total error. Its longitudinal contribution is −0.2819 eV/Å and transverse norm 0.1453 eV/Å, consistent with the −0.3227 eV/Å pair-separation error. Local RMSE is 0.1186 versus 0.0863 outside 4.5 Å. V15/V16 matched reports also repeatedly identify this Ti as the largest local residual on related Ag–Ti registry references.

These are local residual patterns, not proof of a unique chemical cause. Error values can be sensitive to a small number of marked atoms, and the registry references share parent motifs.

## Ranked explanations, evidence, and confounds

| Rank | Candidate factor | Evidence | Counterevidence / confound | Confidence |
|---|---|---|---|---|
| 1 | **Insufficient local response coverage**, especially Ag–Si transverse motion and C/Ti framework neighbors. | Largest persistent residuals are localized to marked Ag/Si and nearby shells; local RMSE is above outside-region RMSE on V16 AgC and AgTi. Adding verified data across V12–V16 coincides with lower matched force error; V16’s six added labels coincide with a drop in all three interface aggregates. | Data, model checkpoint, and epochs jointly changed; only six related local structures are scored. Correlation does not establish which labels caused improvement. | Moderate; strongest actionable hypothesis, not proven cause. |
| 2 | **Inherited label-method/provenance variation.** | V15 removed five uncertain force-only frames, yet V16 still inherits 14 partial-evidence and 12 declaration-only provenance rows; physical method consistency remains unknown. Mixed LCAO/PW-PBE histories are a real risk. | V15→V16 matched forces improve although inherited uncertainty remains; therefore provenance inconsistency alone cannot explain the trend or all residuals. | Moderate as a data-quality risk; low as sole explanation. |
| 3 | **Validation lineage/reuse makes apparent generalization optimistic.** | V14 AgTi test exactly duplicates a training geometry. V15 probes are small perturbations of prior parents. Older membership cannot be reconstructed completely. | This limits evidence about generalization; it does not itself create a force residual on a correctly paired DFT reference. | High as a validation limitation; low as direct error cause. |
| 4 | **Energy/force metric and masking issues may have obscured earlier diagnosis.** | V12 audit showed the training energy RMSE included four force-only rows as zero-energy targets unless masked; correct audit produced a very different energy statistic. Historical revisions also fixed evaluator issues. | Current matched scorer uses explicit energy-per-atom, vector-force, marked-separation, and per-atom errors; V15/V16 published values were hash-checked/reproduced. This does not explain current residuals by itself. | High that old summaries need care; low as current model-error cause. |
| 5 | **Loss balance, optimizer, or architecture.** | The unchanged V12–V16 recipe still misses force gates. A fixed-data controlled ablation could determine whether the force/energy objective or selection contributes. | The force weight is already 1000 vs energy 100; the recipe is stable and force error improves as data expands. No controlled parameter experiment isolates weights, schedule, architecture, or selected epoch. | Low until tested. Do not guess by raising force weight/epochs. |

## V17 plan for A

### Parameters to preserve for the first V17 comparison

Use the V16 baseline unchanged for the first data-only comparison: MACE-MP-0b3 medium, seed 45, CPU, batch size 1, at most 80 epochs, learning rate 1e-4, weight decay 5e-7, energy weight 100, force weight 1000, and estimated E0s. Select only from training/validation. Keep V16 as an immutable baseline. The evidence supports testing targeted data coverage before changing optimizer/loss settings; it does **not** prove that these are globally optimal parameters.

### What the eight paired labels can and cannot answer

The registered V17 ± pairs target (i) Ag–C framework-neighbor response, (ii) Ag–Si transverse Ag response, (iii) Ag–Si framework-Si response, and (iv) Ag–Ti framework-neighbor response. Once both signs for each pair pass archive and atom-ID checks, include those eight labels in the V17 training candidate and compare against V16 with all training controls fixed. Paired signs can show local response sensitivity and whether the sign of the force residual is captured around those four environments. They cannot alone prove generalization to another registry family, separate inherited method effects, or replace a fresh independent test set.

### Minimal discriminating experiment and gates

1. Finish the existing A-positive and B-negative DFT labels independently; verify input/output hashes, IDs, coordinates, convergence, finite values, and force/energy binding before integration. Do not infer any missing sign from symmetry.
2. Freeze the V17 manifest and split after integrating only approved labels. Preserve source method, DFT method, original verification evidence, parent geometry identity, and signed displacement. Enforce lineage-level split isolation; do not use the eight diagnostic structures as a validation or blind set.
3. Train one V17 baseline with the unchanged V16 recipe above. Keep the selected checkpoint blind-label-free. Score the six prior references as a retrospective trend only.
4. Before calling V17 validated, lock a **new independent geometry family** per interface, not a small displacement of a training parent. Score each only after checkpoint selection and report energy, full force-vector RMSE, signed/absolute separation error, top atom vectors, projections, and local/outside errors.
5. Keep the existing gates: ≤10 meV/atom, ≤0.05 eV/Å full vector RMSE per interface, and ≤0.10 eV/Å separation per marked pair. All required gates must pass on independent structures before any MD/TTM stage. A passing separation gate cannot offset a force-vector failure.
6. If V17 still fails, run a **single-factor fixed-data ablation**: hold V17 labels, split, seed, architecture, and schedule fixed and change only one declared force/energy loss-balance setting against the V16 baseline. Report train/valid curves and the same locked DFT references. Choose the factor and candidate value before viewing the new blind scores; do not change weight, epoch limit, and data together.

**Recommendation:** keep the current V16 recipe for V17’s first comparison, finish all eight targeted signed labels, and add a separate locked validation family. If V17 misses a gate, use the resulting matched residuals to decide whether the next experiment should target more coverage or a one-factor loss-balance test. No current evidence supports relaxing thresholds or advancing to MD/TTM.

## Evidence index

- Machine-readable version ledger: `analysis_04/version_ledger.csv`, `analysis_04/version_ledger.json`.
- Same-reference per-structure metrics and identities: `analysis_04/common_geometry_scores.csv`, `analysis_04/common_geometry_scores.json`.
- Driver: `analyze_v9_v16_common.py`.
- Earlier model/version history and V9–V13 score roles: `../../../../periodic_interface_v4/MACE_VERSION_HISTORY.md`, `../../../../WORK_LOG_2026-10.md`, `../v9-v15_force_separation_history.md`, `../v9_v10_v11_force_separation_history.csv`.
- V15/V16 matched per-atom residuals, projections, and metadata repair: `../v16_matched_error_and_metadata_review/report.md` and `../v16_matched_error_and_metadata_review/analysis_04/`.
- V15 inherited provenance limitation: `../v15_provenance_aware_builder_revision/` and `../v15_inherited_method_provenance_trace/` where present in the repository history.
- Gate definitions and staged progression: `../../../../coordination/ITERATION_PROTOCOL.md`.

Earlier failed driver attempts are preserved without being treated as scores: `analysis_01/diagnostics_failure.json` through `analysis_03/diagnostics_failure.json`. The final complete run is `analysis_04/`.
