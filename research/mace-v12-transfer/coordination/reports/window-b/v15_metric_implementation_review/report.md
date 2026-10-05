# V15 metric implementation review — Window B

Reviewed integrated `mace_periodic_v15_interface_energy/evaluate_v15.py` and `preflight.py` after assignment `be21fe5`. Only source code was read. No real reference labels, blind inputs/results/summaries, model/checkpoint files, Torch/MACE/GPAW imports, training or inference were used. Eighteen synthetic mathematical/guard assertions passed. This review does not evaluate the trained model or interrupt A's training.

## Main metric logic

| Item | Finding / formula |
|---|---|
| Energy units | Correct: abs(Epred−Eref)×1000/N gives meV/atom from eV/cell (evaluate lines59). Synthetic 0.48eV/48 atoms =10meV/atom. |
| Force-vector RMSE | Correct: sqrt(mean_atoms(sum_xyz(deltaF²))), in eV/Å (line60). This is sqrt3 times component RMSE for the same array; do not compare thresholds against a component RMSE. |
| Aggregate force RMSE | Correct: sqrt(sum(N_s×RMSE_s²)/sum(N_s)), atom weighted (lines12–17). Synthetic unequal-size example gives sqrt0.0325, not arithmetic mean of structure RMSE. |
| Maximum force error | Correct maximum norm of per-atom error vector, not maximum absolute Cartesian component (line61). |
| Aggregate energy | Mean of per-structure meV/atom errors; intentionally different from a pooled atom-weighted energy statistic. Label/document it as per-structure MAE. |
| Separation sign | Correct cluster convention: `(F_X−F_Ag)·unit(r_X−r_Ag)`; positive means separation. Ag and non-Ag endpoints are chosen by chemical identity, not array order. Reversed synthetic array order gives the same positive separating force (lines51–57). |
| Marked IDs | Preflight requires exact ordered elements, PBC, persistent IDs and central_pair array against the frozen input, unique IDs, coordinate/cell tolerance; code cannot silently remap a wrong ID (preflight32–37). |
| Groups | Historical test and fresh local-registry probes stay separate; interface gates use only fresh probes. Historical scores are regression evidence with acknowledged correlation, not independent transfer certification (evaluate41–42,64–76). |
| Gates | Inclusive ≤10meV/atom, ≤0.05eV/Å vector RMSE, ≤0.10eV/Å separation. Equality accepted; next representable number above threshold rejected; missing metric not accepted. One fresh structure per AgC/AgSi/AgTi is required. |
| Overall status | A failed provisional interface gate gives FAIL. All narrow gates met gives UNDETERMINED, never automatic model PASS; consistent with limited screening scope. |
| Fail guards | Missing fields, nonfinite references/predictions, wrong IDs/PBC/positions and incomplete archives fail. Synthetic extracted guards demonstrate reference and geometry rejection without any real data reads. |

## Definite robustness defects and proposed patch

### D1 — derived nonfinite values fail late (low severity)

`evaluate_v15.py:47` checks finite model energy and force components, but lines56–61 form differences, squared norms and projections without checking the derived values. Finite synthetic deltaF=[1e308,0,0] overflows its squared norm to inf. The later `json.dumps(..., allow_nan=False)` at line78 rejects inf, so this does **not** publish a valid-looking metric PASS. However, it fails late, after processing later frames and creating the output directory, without the offending-frame metric diagnostic. A similar subtraction overflow is possible from finite extreme energies. The report's synthetic case reproduces the exact norm expression and JSON rejection; no real predictions were generated.

Proposed fix: check all derived per-frame metrics for finiteness immediately after constructing the row, with its label in the error; check aggregate results as well, since weighted accumulation can overflow. In ordinary physical ranges this is a robustness change; this review does not claim current V15 predictions overflow.

### D2 — selection record hash can describe a different snapshot (medium evidence-integrity severity)

Line27 parses the selection record once. Line73 later hashes its path again after potentially long inference. If the file is rewritten between those operations, selection checks use the original object but output records the SHA of replacement bytes. The final model SHA check (line72) protects model bytes, not the selection-record snapshot. The counterexample is a reviewed initial record for selected-A followed by a replacement for selected-B or a changed blind-selection declaration; the output SHA can describe the unvalidated replacement. This is a source-level two-read race; it is not evidence that A actually rewrote its record.

Proposed fix: read selection bytes once, calculate/store their SHA, parse those bytes, and require the file still matches that SHA before publishing results. Report the checked initial SHA. A should keep completion evidence immutable as well; the current design treats A's signed-off declarations as review evidence rather than parsing an epoch log automatically.

`proposed_patch.diff` addresses D1/D2 only and passed static compilation as modified source in memory. It was not applied to production or executed with a model. The compact patch uses zero context; A may review it manually or apply it with the appropriate zero-context patch option after review. A reviews/integrates it. This task is a code review; finding defects is a completed review outcome, not a failed scientific calculation.

## Limitations and policy decisions (separate from defects)

1. **Periodic-image convention:** lines55–57 explicitly retain V14's direct cluster coordinate vector. Preflight requires PBC identity but does not enforce a minimum-image pair vector or verify that the marked pair lies in a common chosen image. Synthetic periodic positions Ag=(0.1,0.1,0), X=(9.9,1.1,0) in a10Å cell give direct displacement(9.8,1,0), MIC displacement(−0.2,1,0). With force difference(1,0,0), projections differ (values in checks.json), and even absolute projected errors can differ. This is a documented cluster-convention limitation for wrapped/extended-cell reuse, not proof any frozen V15 geometry or actual score is wrong. Do not silently change the V14/V15 comparison definition to MIC; A should declare and verify the intended image convention before extending scope.
2. **Selection evidence:** a reviewed selection record must declare training_exit_code0, completed_epochs80, selected_epoch in1..80, training/validation-only selection and no blind selection, with model and completion evidence hashes. This intentionally trusts A's manual completion review. Hashes alone do not prove80 epochs occurred. A should state whether selected_epoch is1-based or translates zero-based log epochs0..79; an unconverted zero-based epoch0 would be rejected. This is a record-contract issue, not evidence the current selection is wrong.
3. **Archive evidence:** static preflight code requires all three archives PASS, declared checks true, input linkage, complete/converged state, four ranks, archive member hashes and source/output geometry before returning references. No archive runtime branch was executed. It trusts the committed verifier's PASS/hash chain, rather than independently proving original SCF physics at review time.
4. **Sample scope:** each fresh interface currently has one probe; energy MAE and RMSE across one probe are narrow tests. No morphology, thermal/liquid or broad periodic-interface transfer follows. Fourteen partial and twelve explicit inherited training declarations still lack full original-run revalidation.
5. **Exact boundaries:** floating-point arithmetic near a decimal gate may round above it. The implementation uses exact inclusive comparison without extra tolerance; do not introduce an undocumented tolerance based on model outcomes.

## Reproduction and publication

From repository root:

```bash
/workspace/.venvs/gpaw-mpi/bin/python -B research/mace-v12-transfer/coordination/reports/window-b/v15_metric_implementation_review/synthetic_checks.py
sha256sum -c research/mace-v12-transfer/coordination/reports/window-b/v15_metric_implementation_review/SHA256SUMS.txt
```

The check script extracts only aggregate/require/labels/geometry AST functions, runs fabricated arrays/objects, and hashes the two reviewed source files. It never imports or calls evaluation main/preflight checks. `checks.json` records exact reviewed hashes,18 assertions, two robustness findings and the periodic example. Historical/current scientific result files are not dependencies. Production files and A's training/owner records remain unchanged.
