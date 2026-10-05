# V14/V15 force and separating-force driver diagnosis

## Scope and status

This is a matched retrospective comparison: both frozen models were evaluated on the same three already-scored V15 registry probes. CPU float64 inference used V14 model `41bd411355b034b0dab85dbdbda932ecfc394122cb67a2b51a550bd04207645d` and V15 model `f258bee17bc46e9399c1f1d13bb5ed30d4a03997b33baac53b5e58c2df6ebc1a`. The exact V15 scores reproduce the frozen evaluation within 1e-06 eV/Å or eV/cell. All probe input/reference hashes, IDs, formula, coordinates, cell, PBC, and SCF verification were checked. The three V15 geometries are no longer blind; their labels are acquisition/history. No V16 result/reference, new DFT, training, MD, or TTM was used.

## Common geometry comparison

| Interface | V14 force vector RMSE | V15 force vector RMSE | V15−V14 | V14 separation error | V15 separation error | V15−V14 |
|---|---:|---:|---:|---:|---:|---:|
| AgC | 0.142882 | 0.135282 | -0.007600 | 0.315007 | 0.387440 | +0.072433 |
| AgSi | 0.335893 | 0.313283 | -0.022611 | 0.012284 | 0.036305 | +0.024021 |
| AgTi | 0.116872 | 0.117514 | +0.000642 | 0.307976 | 0.350126 | +0.042150 |

Energy MAE remains small on these geometries; exact model-by-model energy, component RMSE and max-vector scores are in `common_geometry_scores.csv`. Deltas are descriptive on a shared input, not causal estimates. The probes share small-cluster source motifs that appear in both training sets; this is a useful controlled geometry comparison but weak evidence for transfer to a different morphology.

## What changes on the same structures

- **AgC:** force RMSE 0.1429 → 0.1353 eV/Å; separation absolute error 0.3150 → 0.3874 eV/Å.
- **AgSi:** force RMSE 0.3359 → 0.3133 eV/Å; separation absolute error 0.0123 → 0.0363 eV/Å.
- **AgTi:** force RMSE 0.1169 → 0.1175 eV/Å; separation absolute error 0.3080 → 0.3501 eV/Å.

V14→V15 train edits were five removals and three additions: four force-only LCAO rows and one unresolved AgTi row left; three directly archived residual-shell labels entered. V15 train size is 29 versus V14's 31. Validation and test files are byte-identical. V15 preserved inherited method declarations where available, but the provenance review found zero fully confirmed inherited original calculations among its 15 recovery targets; 14 are partial and one unresolved row was excluded. V15's residual references have direct runner/log/archive evidence. These dataset/provenance edits happened together, so the matched delta cannot identify which edit helped or hurt.

Both training runs used the same declared CLI settings: seed 45, CPU, batch size 1, 80 epochs, learning rate 1e-4, weight decay 5e-7, energy weight 100 and force weight 1000, with the same foundation-model SHA. The selected models also have matching serialized state-dictionary key/shape schemas. Checkpoints were selected from training/validation, V14 one-based epoch 71 and V15 one-based epoch 75. V15's selected validation force RMSE was 0.03882 eV/Å; the V14 selected value was 0.04410. At epoch 79 those values were 0.04535 and 0.04073. Validation improved slightly, while this matched local registry set still fails force gates; validation does not represent every local transverse response.

## Spatial and projection diagnosis

Forces are partitioned without overlap into the marked Ag/X pair, other Ag atoms within 4.5 Å MIC of either marked atom, non-Ag framework atoms in that same shell, and all atoms outside it. This 4.5 Å cutoff matches the previous localization analysis; it was not selected by optimizing these results. Each category includes vector and x/y/z component RMSE, species breakdown and SSE share in `force_decomposition.json`; top-error atoms include persistent IDs and nearest MIC neighbors.

### AgC (Ag#688433–C#604804)

| Model | Pair Ag signed term | Pair X signed term | Sum separation error | Ag transverse norm | X transverse norm | Local Ag-neighbor RMSE | Local framework RMSE | Outside RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| v14 | -0.0551 | -0.2599 | -0.3150 | 0.3871 | 0.0756 | 0.0542 | 0.1506 | 0.0877 |
| v15 | -0.0816 | -0.3058 | -0.3874 | 0.3751 | 0.1062 | 0.0780 | 0.1358 | 0.0825 |

V15 top five atom errors: Ag#688433 0.384 eV/Å; C#604804 0.324 eV/Å; Ti#602709 0.315 eV/Å; C#604602 0.270 eV/Å; Ti#604589 0.241 eV/Å.

### AgSi (Ag#683609–Si#366934)

| Model | Pair Ag signed term | Pair X signed term | Sum separation error | Ag transverse norm | X transverse norm | Local Ag-neighbor RMSE | Local framework RMSE | Outside RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| v14 | +0.0288 | -0.0410 | -0.0123 | 0.1694 | 1.2640 | 0.5502 | 0.2198 | 0.1904 |
| v15 | -0.0004 | -0.0359 | -0.0363 | 0.1523 | 1.0977 | 0.5189 | 0.2099 | 0.1898 |

V15 top five atom errors: Ag#685172 1.325 eV/Å; Si#366934 1.098 eV/Å; Ag#730721 0.544 eV/Å; C#368838 0.483 eV/Å; C#368645 0.436 eV/Å.

### AgTi (Ag#687247–Ti#479027)

| Model | Pair Ag signed term | Pair X signed term | Sum separation error | Ag transverse norm | X transverse norm | Local Ag-neighbor RMSE | Local framework RMSE | Outside RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| v14 | -0.0557 | -0.2522 | -0.3080 | 0.0329 | 0.1052 | 0.0848 | 0.1506 | 0.1045 |
| v15 | -0.0833 | -0.2668 | -0.3501 | 0.0409 | 0.1217 | 0.0928 | 0.1452 | 0.1021 |

V15 top five atom errors: Ti#479027 0.293 eV/Å; Ag#682693 0.287 eV/Å; Ag#701905 0.284 eV/Å; C#480975 0.222 eV/Å; Ti#479021 0.209 eV/Å.


The marked Ag and X terms sum exactly to the scalar separation error. Their transverse residuals are not included in that scalar. A small radial error can therefore coexist with a large vector error; changing only distance or radial-force weighting would not target that failure. Conversely, a large scalar separation error can have both pair terms aligned in the same direction. Pair distance was checked under both direct cluster coordinates and MIC; those agree for all three probes. Local neighbor distances use MIC, while the legacy signed separation score uses the direct Ag→X vector.

## Ranked hypotheses and decisions for A

| Rank | Evidence and counterevidence | Confidence | Confound | Discriminating test | Suggested action |
|---|---|---|---|---|---|
| 1. Registry/transverse coverage and local environment | Same-geometry decomposition identifies local Ag/framework errors, including transverse components; AgSi's scalar separation can look small while its vector error remains large. The three probes share inherited cluster motifs and only sample a few registry shifts. | Medium-high that this is a concrete error mode; low on its share of the V15→V14 change. | V15 simultaneously changed labels and exclusions; only one geometry per interface here. | On fixed parent structures, compare paired ± in-plane shifts at constant Ag–X normal distance; separately perturb framework coordinates with Ag held fixed. Keep identical DFT method and training recipe. | V16's three registry inputs probe lateral registry, while its three distance inputs probe normal response. They do not test independent framework motion or new morphology. A should review all six labels before deciding inclusion. |
| 2. Historical reference-method uncertainty | Provenance review documents 14 partially recovered inherited labels and one unresolved AgTi row excluded; inherited energy conventions and convergence are not uniformly proven. The three current probes themselves have passing verification and direct matched GPAW archives. | Medium that this is a dataset-risk factor; low that it explains these specific force residuals. | The evaluated probes were newly calculated with the same documented protocol; no paired re-label data exists. | Recover source inputs/logs, or make a separately authorized exact-geometry matched-method reference set with full getter/energy convention and force-order metadata. | Keep current labels unchanged; do not retrofit free energy or infer missing setup. Preserve per-row evidence tiers and isolate convention comparisons. |
| 3. Training objective/optimization or checkpoint generalization | V15 validation force RMSE is modestly better than V14 but held-out registry force errors do not meet the gate. Both runs used same declared settings; no evidence singles out loss weights, learning rate, epoch count or architecture. | Low to medium. | Changed data composition, inherited target quality and validation coverage; single seed and tiny validation split. | After a fixed data/reference baseline is frozen, run one-factor matched ablations (first force-weight 1000 vs 2000 with seed/settings held fixed; select only on a strengthened validation set). | For A's immediate V16 cycle, preserve the current CLI settings so the six new structures' data effect can be read. Do not increase epochs or force weight as an untested repair. |
| 4. Force/vector scoring, atom order or PBC bug | The frozen V15 scores are reproduced; model/reference structures match atom IDs, order, positions, cell and PBC. Vector RMSE and Cartesian component RMSE are separately reported. Pair separation uses the documented direct vector, which agrees with MIC distance for these probes. | Low for these three comparisons; checks rule out common implementation/order errors. | This audit validates these exact archives and scorer convention, not every inherited dataset serialization or every future structure. | Keep automated exact-ID/order and vector-vs-component regression checks in future evaluation. | Retain the existing metric definition; do not replace it with a component RMSE or change the PBC convention silently. |
| 5. DFT energy/smearing/k-point mismatch as direct force cause | Current probes use a documented PBE plane-wave/Gamma/Fermi-Dirac 0.1 eV driver; native extrapolated energy and free energy are separately stored. Historical labels show multiple cells/sampling conventions and incomplete original evidence. Energy errors on current probes are small. | Low for a direct cause on these references; broader target consistency remains unresolved. | Energy offsets do not imply force errors, and no matched convention perturbation exists. | Only if provenance recovery warrants it, compare native/free energy and forces on exact same geometries under documented converged conditions. | Keep REF_energy convention fixed and explicit; record both energies, smearing, k-points, forces, software and convergence for future labels. |

The parameter recommendation is intentionally conservative and directly usable: for the next data-integration cycle, retain seed 45, CPU, batch 1, 80 epochs, lr 1e-4, weight decay 5e-7, energy weight 100 and force weight 1000 while adding only reviewed V16 labels. That keeps the training recipe controlled so A can attribute differences more cleanly to the changed training set. If errors remain after source and split audits, a future fixed-data two-arm force-weight comparison (1000 vs 2000) is a proposal, not an instruction to train now; use additional validation geometries for checkpoint selection and keep a new independent blind set untouched.

## Six planned V16 inputs: what they test

- `AgC_widespan_distance_acq_v16_01` (AgC_wide_distance_p0p35, radial/distance, direct gap 2.610 Å; parent 2.260 Å): Ag-X normal distance response; it cannot resolve independent framework morphology/thermal response; provenance of inherited labels; optimizer causality. Input hash `75841f91e522499e24953b3d6c9454cccc09bab484e1f55001ad26380686dd21`.
- `AgSi_widespan_distance_acq_v16_01` (AgSi_wide_distance_p0p50, radial/distance, direct gap 2.940 Å; parent 2.440 Å): Ag-X normal distance response; it cannot resolve independent framework morphology/thermal response; provenance of inherited labels; optimizer causality. Input hash `cccdfcdd8891216a4078d11cd618dd862daae2c368bfe74cef7c3a73ac3b6049`.
- `AgTi_widespan_distance_acq_v16_01` (AgTi_wide_distance_p0p35, radial/distance, direct gap 2.938 Å; parent 2.588 Å): Ag-X normal distance response; it cannot resolve independent framework morphology/thermal response; provenance of inherited labels; optimizer causality. Input hash `d9a2a87496d133e22778f0531c2a2742a61d0c5a86c61a04fea44d4d500f6115`.
- `AgC_widespan_registry_acq_v16_01` (AgC_wide_registry_p0p45, lateral registry, direct gap 2.304 Å; parent 2.260 Å): lateral Ag registry/contact response; it cannot resolve independent framework morphology/thermal response; provenance of inherited labels; optimizer causality. Input hash `061e1c3374962e2da4ee1378f36f00895df8d669758334620277e275821c051a`.
- `AgSi_widespan_registry_acq_v16_01` (AgSi_wide_registry_m0p30, lateral registry, direct gap 2.458 Å; parent 2.440 Å): lateral Ag registry/contact response; it cannot resolve independent framework morphology/thermal response; provenance of inherited labels; optimizer causality. Input hash `92dd6a42087553d506e438cfae153f2bc58c1bc415d70f3f58da33d0e34cea91`.
- `AgTi_widespan_registry_acq_v16_01` (AgTi_wide_registry_p0p45, lateral registry, direct gap 2.627 Å; parent 2.588 Å): lateral Ag registry/contact response; it cannot resolve independent framework morphology/thermal response; provenance of inherited labels; optimizer causality. Input hash `6c619d9fcda72876449c2792b4d6c4d3db4afcde86bf976cb824e5e15ef74a6a`.

The six are training acquisitions, not blind validation. The parallel-set distances are reported in `audit_evidence.json` against the parent pair value; the main-set structures primarily alter registry. The minimal follow-up proposal is paired ±0.15 Å transverse translations (two independent in-plane axes, fixed normal separation) and paired ±0.04 Å local framework displacements for selected C/Ti/Si neighbors, with the Ag sublattice fixed for the framework arm. Generate matching controls from the same parent and use one unchanged PBE/PW/smearing/k-point/convergence recipe. Assign labels only after A approves. Keep any new independent-morphology/thermal validation inputs separate from these acquisition variants.

## Reproducibility and limits

Run from the repository root with the B MACE CPU environment:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /workspace/.venvs/mace-v14-assist/bin/python -B \
  research/mace-v12-transfer/coordination/reports/window-b/force_separation_driver_diagnosis/analyze_matched_v14_v15.py \
  --repo-root /workspace/- \
  --output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/force_separation_driver_diagnosis/analysis_rerun_01
```

Outputs: `common_geometry_scores.csv`, `force_decomposition.json`, `audit_evidence.json`, this report and `SHA256SUMS.txt`. Evidence hashes cover both selected models, selected/evaluation records, three scored references and verifications, source input manifests/frames, V14/V15 split and training logs/scripts, inherited provenance/convention audits, V16 input manifests/frames, and the exact runner. V16 calculation directories are not read. No new data labels, model fitting, or model selection were produced.
