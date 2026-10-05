# V15 per-atom force-error localization

Selected/scored model SHA256: `f258bee17bc46e9399c1f1d13bb5ed30d4a03997b33baac53b5e58c2df6ebc1a`. Exactly three already-scored fresh V15 probes were analyzed using CPU float64. Post-evaluation authorization in B_CONTINUOUS_WORK.md permits these references and inference; no other unseen labels/models were inspected. No model selection, tuning, training, DFT or MD/TTM occurred. V15 remains **FAIL**.

## Reproduced scores and spatial distribution

Local = within4.5Å of either marked atom using minimum-image distances. Separation uses the direct Ag→X cluster vector, matching A’s scorer. Cutoff is the historical localization cutoff, not tuned to V15 outcomes. Species/Cartesian details and every atom’s ID/error are in JSON/CSV.

| Interface | All vector RMSE | Local RMSE | Outside RMSE | Local error SSE share | Maximum atom error | Separation error |
|---|---:|---:|---:|---:|---:|---:|
| AgC | 0.135282 | 0.158746 | 0.082481 | 86.1% | 0.383842 | 0.387440 |
| AgSi | 0.313283 | 0.386218 | 0.189834 | 83.4% | 1.325495 | 0.036305 |
| AgTi | 0.117514 | 0.128688 | 0.102142 | 66.1% | 0.293258 | 0.350126 |

Units are eV/Å except SSE share. All per-probe vector RMSE/max, separation error and predicted energy reproduce published values within1e-6eV/Å or1e-6eV/cell. Pooled fresh-probe metrics also reproduce A: force-vector RMSE 0.210542577eV/Å; energy MAE 2.269230meV/atom; separation max 0.387439909eV/Å. Exact numerical deltas are in JSON. This tolerance accommodates CPU float64 library/order rounding and is not a relaxation of any gate.

## Contact projections versus full-vector errors

| Interface | Ag ID | Contact ID | Ag signed separation-error contribution | X signed contribution | X transverse error norm |
|---|---:|---:|---:|---:|---:|
| AgC | 688433 | 604804 | -0.081602 | -0.305838 | 0.106217 |
| AgSi | 683609 | 366934 | -0.000409 | -0.035896 | 1.097659 |
| AgTi | 687247 | 479027 | -0.083300 | -0.266826 | 0.121672 |

The two signed contributions sum to predicted-minus-reference separation error. A small scalar projection does not constrain transverse forces or errors on other atoms.

### Ag–Si: small separation error masks large transverse/local errors

The marked Si#366934 has full-vector error1.098246eV/Å, of which1.097659eV/Å is transverse to the marked separation direction. Its longitudinal contribution is only−0.035896eV/Å; Ag#683609 contributes−0.000409eV/Å, giving total−0.036305eV/Å. The small scalar error is chiefly a projection issue here, not large longitudinal cancellation. Ag#685172 is the largest-error atom at1.325495eV/Å, not the marked Ag. Ag and Si contribute47.9% and29.1% of force-error SSE respectively. Local RMSE0.386218 and outside RMSE0.189834 both exceed the0.05 gate. Broad Ag/framework transverse coverage is needed; fixing only the marked scalar distance is insufficient.

### Ag–C

86.1% of force-error SSE lies in the local shell. Highest errors include Ag#688433(0.383842), C#604804(0.323758), Ti#602709(0.315083), C#604602(0.269516) and Ti#604589(0.241081)eV/Å. Contact C dominates the longitudinal separation error(−0.305838), while marked Ag contributes−0.081602 and also has transverse error0.375068. Carbon/Ti shell errors and lateral Ag placement both matter.

### Ag–Ti

66.1% of force-error SSE is local; the remaining33.9% outside-shell share is non-negligible. Marked Ti#479027 is the largest error at0.293258eV/Å. Its longitudinal contribution−0.266826 and Ag#687247 contribution−0.083300 give−0.350126eV/Å total separation error. The model predicts negative separation(−0.122405) where DFT is positive(+0.227722). Other large errors include Ag#682693, Ag#701905, C#480975 and Ti#479021. Distance/contact coverage and outside framework responses both need attention.

## Species breakdown

| Interface | Element | Atoms | Vector RMSE eV/Å | Error SSE share |
|---|---|---:|---:|---:|
| AgC | Ag | 2 | 0.276961 | 17.5% |
| AgC | C | 13 | 0.148897 | 32.8% |
| AgC | Si | 10 | 0.036408 | 1.5% |
| AgC | Ti | 23 | 0.135707 | 48.2% |
| AgSi | Ag | 17 | 0.375555 | 47.9% |
| AgSi | C | 13 | 0.237891 | 14.7% |
| AgSi | Si | 8 | 0.426579 | 29.1% |
| AgSi | Ti | 13 | 0.178946 | 8.3% |
| AgTi | Ag | 28 | 0.095685 | 37.9% |
| AgTi | C | 11 | 0.126195 | 25.9% |
| AgTi | Si | 1 | 0.076637 | 0.9% |
| AgTi | Ti | 9 | 0.163047 | 35.4% |

## V14 pattern context — different probes, not matched model comparison

| Interface | V14 local/outside RMSE | V15 local/outside RMSE |
|---|---|---|
| AgC | 0.1762 / 0.0937 | 0.1587 / 0.0825 |
| AgSi | 0.1266 / 0.1079 | 0.3862 / 0.1898 |
| AgTi | 0.1775 / 0.0834 | 0.1287 / 0.1021 |

Both historical and current patterns show local errors and nonzero outside-shell errors. These are different registry probes and differently trained models: the table must not be read as a controlled V14→V15 improvement/degradation measurement. A matched comparison would need an explicitly designed common evaluation set; none was executed here.

## Existing candidates: useful axes and remaining gaps

- **AgSi_wide_registry_m0p30** (backup existing passing +0.20 registry): prioritize A review of lateral Ag/environment response because Si/Ag transverse errors dominate. The rigid shift probes one lateral axis; it does not establish coverage of arbitrary Si displacement, independent morphology or every high-error Ag atom. **AgSi_wide_distance_p0p50** extends the normal distance bracket but does not directly target the transverse Si error. Keep it as complementary coverage, not the sole remedy.
- **AgC_wide_registry_p0p45** probes lateral Ag-shell coupling; **AgC_wide_distance_p0p35** broadens contact separation. Neither moves C/Ti framework atoms independently, so the large C/Ti shell errors may require a later controlled framework-acquisition design after A review.
- **AgTi_wide_distance_p0p35** broadens normal contact response and **AgTi_wide_registry_p0p45** samples lateral coupling. Neither alone addresses the33.9% outside-shell SSE or independent Ti/C framework motion.
- Existing ±0.12 surrounding-shell proposals remain near-parent correlated variants and were not among the six PASS/distinct initial selections. Their original perturbation IDs do not necessarily match the newly largest-error IDs; do not call them targeted repairs or promote previously flagged contacts based on this report. No candidate geometry, priority manifest, owner assignment or launch block was changed.

The existing six-candidate pack remains blocked pending A review and explicit ownership/production integration. The observations guide future acquisition/history, not further tuning against an unseen validation set. Already-scored V15 probes are no longer blind; if A later uses their labels in training, their future role is acquisition/history, and a new independently withheld set is required for unseen validation. This report does not create new blind tests or authorize DFT, retraining, MD/TTM.

## Artifacts and reproduction

`analysis_01/localization.json`: per-interface/species/region metrics, signed projections, top10 IDs with nearest4 MIC neighbors, exact score deltas, runtime versions and model/reference/evidence hashes. `analysis_01/per_atom_errors.csv`: all148 atoms with Cartesian DFT/predicted/residual forces, vector error, IDs and local classification.

The existing `/workspace/.venvs/mace-v14-assist` provides Torch2.5.1+cpu, mace-torch0.3.16 and ASE3.29.0; no installation or environment configuration change was needed. One CPU thread was used.

Actual executed command:

```bash
cd /workspace/-
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /workspace/.venvs/mace-v14-assist/bin/python -B \
  research/mace-v12-transfer/coordination/reports/window-b/v15_force_error_localization/localize_v15.py \
  --repo-root /workspace/- \
  --output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v15_force_error_localization/analysis_01
```

The output is now nonempty; repeat execution refuses overwrite. Use a new empty report child only if a rerun is authorized/necessary. On inference/hash/metric failure, preserve diagnostics_failure.json and stop; do not substitute another checkpoint/reference or weaken tolerances.

From repository root verify `SHA256SUMS.txt`. It covers compact deliverables and accessed evidence, including a hash pointer to the existing selected model; no large model/checkpoint is copied to B reports.
