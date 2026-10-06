# Matched V15/V16 force-error review

## Scope and checks

Both window-A-reviewed frozen checkpoints were run on the same six archived, already-scored registry references: three originally scored with V15 and three originally scored with V16. This is retrospective diagnosis, not fresh validation or model selection. No DFT, training, MD or TTM was run.

The selected checkpoints were V15 epoch 75 (`f258bee17bc46e9399c1f1d13bb5ed30d4a03997b33baac53b5e58c2df6ebc1a`) and V16 epoch 80 (`918aad24341c6892d65ec0d31b19baf610ac0efde1cfd9d2e1c0567320e0cf98`). Both training entry scripts specify the same foundation model, seed 45, CPU, batch size 1, 80 maximum epochs, learning rate 1e-4, weight decay 5e-7, energy weight 100 and force weight 1000. V16 adds three wide-span registry and three wide-span distance labels. The training recipe is controlled, while the selected epoch and six training labels differ together; this comparison therefore supports a data/fit association, not a causal estimate for the six labels alone.

Before inference, the V16 metadata repair was rechecked: both provenance totals equal 35; train/valid/test are 35/2/3; current pins match; the stored pre-repair manifest and pins match commit `b6ac255`; and split file bytes are unchanged from that commit. V16 preflight passed. The six reference hashes, archive verification records, input hashes, atom identities, PBC, cell and coordinates matched the original evaluation records. Each model reproduced its own three published scores within 1e-6 force units and 1e-6 eV/cell. Full inputs and score-reproduction deltas are in `analysis_04/matched_comparison.json`.

The local shell uses ASE MIC distances to either marked Ag/X atom with the fixed 4.5 Å radius. The separation direction follows the original scorer's direct Ag-to-X vector. Inference used CPU float64 and one Torch thread. “Framework” means every non-Ag atom.

## Same-geometry results

All six force RMSE values decreased with V16 on the exact same reference geometry. The pooled vector RMSE over 296 atoms decreased from 0.1852 to 0.1464 eV/Å (21%). This is a consistent matched improvement, but every interface-level V16 score on its own fresh holdout still exceeds the 0.05 eV/Å gate.

| Reference geometry | Role | V15 force RMSE | V16 force RMSE | Change | V15 separation error | V16 separation error |
|---|---|---:|---:|---:|---:|---:|
| AgC V15 holdout | V15 | 0.1353 | 0.0722 | -0.0631 | 0.3874 | 0.1027 |
| AgSi V15 holdout | V15 | 0.3133 | 0.2804 | -0.0329 | 0.0363 | 0.0079 |
| AgTi V15 holdout | V15 | 0.1175 | 0.0729 | -0.0446 | 0.3501 | 0.1583 |
| AgC V16 holdout | V16 | 0.1831 | 0.1295 | -0.0536 | 0.0285 | 0.0417 |
| AgSi V16 holdout | V16 | 0.1334 | 0.0961 | -0.0372 | 0.0260 | 0.0793 |
| AgTi V16 holdout | V16 | 0.1478 | 0.1053 | -0.0426 | 0.6236 | 0.3227 |

All force values and changes are eV/Å. Separation error is absolute eV/Å. Per-interface pooled force RMSE across each interface's two geometries:

| Interface | V15 | V16 | Relative change |
|---|---:|---:|---:|
| Ag–C | 0.1610 | 0.1048 | -34.9% |
| Ag–Si | 0.2408 | 0.2096 | -12.9% |
| Ag–Ti | 0.1335 | 0.0905 | -32.2% |

Ag–Si remains the largest force-error source. The V16 model still has a 0.1053 eV/Å force RMSE on its Ag–Ti holdout, and its Ag–Ti separation error there is 0.3227 eV/Å. A pooled retrospective score cannot override these per-interface gates.

## Where the residual errors sit

- **Ag–Si transverse response persists.** On the same AgSi V15 geometry, V16 improves full RMSE only from 0.3133 to 0.2804. The marked Si atom's transverse error grows from 1.098 to 1.135 eV/Å, while the scalar separation error shrinks from 0.0363 to 0.0079. The latter results from a small net longitudinal projection; it does not mean the contact force vector is accurate. Ag#685172 is the maximum-error atom in both models (1.325 to 1.344 eV/Å). On the V16 AgSi geometry, V16's worst atoms are Ag#710704, Ag#685172 and Si#367120; errors are spread across the local shell and outside atoms such as C#368645. This points to transverse Ag/framework response coverage, not only the marked-pair distance.
- **Ag–Ti has a longitudinal contact failure.** On the V16 AgTi geometry, V16 predicts -1.3498 eV/Å against -1.0271 eV/Å reference separation force. The marked Ti error is 0.3172 eV/Å, including a -0.2819 eV/Å longitudinal contribution; marked Ag adds -0.0407 eV/Å. This is not primarily a transverse projection artifact. The V16 matched local/outside RMSE is 0.1186/0.0863 eV/Å, so framework response beyond 4.5 Å also remains material.
- **Ag–C scalar success masks lateral Ag error.** On the V16 AgC geometry, V16 passes the separation gate (0.0417) but the marked Ag#688433 vector error is 0.4451 eV/Å, of which 0.4305 is transverse. Local RMSE remains 0.1507 versus 0.0830 outside. A radial metric alone would miss this.
- Across all six references, V16 improves both local and outside-shell RMSE on every matched geometry. The residual distribution remains mostly local (about 62–93% of force-error SSE per geometry), but outside-shell errors are not negligible, particularly for AgSi and AgTi.

## Evidence and next diagnostic

The matched six-geometry comparison shows V16's selected model responds better than V15 across the tested fixed geometries, with the largest interface-level gains in Ag–C and Ag–Ti. V16's added labels and the same training recipe are consistent with improved coverage; because the data and selected epoch changed together, the comparison does not isolate either the six labels or loss weights as the cause. The small six-geometry set also cannot establish transfer to a new morphology or thermal distribution.

Use the already registered paired transverse/framework proposal pack for the next controlled acquisition review. It separates AgSi transverse Ag movement from Si-contact movement, and C/Ti contact response from lateral Ag registry. The currently registered V17 opposite-sign pairs are the compact first diagnostic set; do not expand the label list before comparing those results. The most informative priorities are:

1. AgSi transverse Ag and framework-Si sign pairs, because its transverse marked-Si error stays large despite a passing radial error.
2. AgTi framework-Ti sign pair, because the V16 holdout shows a large longitudinal Ti-contact miss plus nonzero outside-shell error.
3. AgC framework-C sign pair, while tracking the marked Ag's transverse vector error separately from separation force.

To separate coverage from loss weighting, first compare a data-only retrain with unchanged training settings against V16 on the same fixed references. Then, only if needed, run a one-factor loss-weight ablation on the same fixed dataset and seed schedule. Both comparisons need a newly frozen independent reference set; these six already-scored probes are no longer fresh blind validation. A's active diagnostics continue independently. B will complete the queued V9–V16 retrospective before starting its queued V17 acquisitions, per the updated priority instruction.

## Artifacts

- `analysis_04/matched_comparison.json`: hashes, metadata checks, every geometry/model result, local/species regions, marked-pair projections, top atoms with MIC neighbors, aggregate comparison and proposal mapping.
- `analysis_04/matched_metrics.csv`: compact six-row paired score table.
- `analysis_04/per_atom_errors.csv`: all atoms for both models on all six geometries.
- `analysis_04/recommendations.json`: next acquisition and controlled-comparison recommendations.
- `analyze_matched_v15_v16.py`: guarded reproduction script. Earlier failed diagnostic attempts are preserved in `analysis_01`–`analysis_03`; the verified output is `analysis_04`.
