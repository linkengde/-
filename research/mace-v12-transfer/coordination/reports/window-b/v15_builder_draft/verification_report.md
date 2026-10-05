# V15 builder draft exercise and provenance report

Audited source snapshot: main `91287a415eba8d3c69dbef53b64cf5fb2c9057de`. Exercised script SHA256: `f9a1d8c81cb80a006f87610be04abf594fe616332b5517d38400a8b56023fb2a`. Generated data is confined to this report directory; production data/builders, A tasks, frozen inputs, models, training settings, worklog and global manifest were not changed.

**Engineering checks PASS; physical-method consistency and V15 validation remain UNKNOWN.** A integration/review is pending. The draft is not approval to train. No V15 holdout labels, summaries, DFT outputs or active A files were read.

## Actual generated data

| Split | frames | energy labels | force labels | disposition |
|---|---:|---:|---:|---|
| train | 30 | 30 | 30 | 31 − 4 + 3 = 30; finite labels validated |
| valid | 2 | 2 | 2 | V14 bytes copied verbatim; historical evaluation only |
| test | 3 | 3 | 3 | V14 bytes copied verbatim; historical evaluation only |

The energy composition matrix in `[Ag,C,Si,Ti]` order has rank **4/4**; matrix rows and singular values are in the manifest. All frame labels are read from their existing REF fields; no absent energies are filled with zeros. Counts are checked before and after serialization. Atom order, elements, persistent IDs where present, positions, cell/PBC and energy/force values survive output verification. Historical config tags can be stale: a `force_only` tag is not used to infer a missing energy when a committed REF energy exists; its method provenance still remains unknown.

## Explicit exclusions and verified additions

| Excluded V14 frame | config | reason |
|---|---|---|
| 9 | `AgC_primary_2.0A_LCAO` | force-only LCAO with missing original source/basis/boundary/convergence provenance |
| 10 | `AgC_primary_2.4A_LCAO` | force-only LCAO with missing original source/basis/boundary/convergence provenance |
| 11 | `AgTi_2p1A_LCAO` | force-only LCAO with missing original source/basis/boundary/convergence provenance |
| 12 | `AgTi_2p7A_LCAO` | force-only LCAO with missing original source/basis/boundary/convergence provenance |

Exactly the three AgC/AgSi/AgTi residual-shell acquisitions enter training. For each, the draft checks input and source-parent hashes against the acquisition manifest; archive PASS and every declared member hash; summary/progress convergence; four MPI ranks; PW-PBE/500 eV/Gamma metadata; finite energy and `(N,3)` forces; summary energy/counts; and input/output ID order, elements, positions, cell/PBC and marked-pair tags. An exact duplicate of an already retained training geometry aborts. All checks passed. The 33 wide-span proposals are not loaded or included.

## Frozen V15 input isolation and historical leakage

Only the V15 holdout input manifest and its three explicit input files are permitted by the read guard. Input hashes and absence of labels are checked. Exact/near geometric overlap against train, valid and test is **zero**. The guard blocks holdout calculation/output paths before content is read. The frozen input geometry and all its labels remain outside training/validation/checkpoint/hyperparameter selection.

Historical AgTi `AgTi_2p7A_LCAO` train and `AgTi_2p70_periodic_PW_PBE_energy_force_holdout` test are identical by IDs/elements/cell/centered coordinates (RMS/max 0), but PBC differs. Complete calculator-input identity is unknown and is not inferred from geometry. The LCAO train row is excluded; the test bytes are retained and explicitly marked non-independent historical regression evidence.

After exclusion, three retained train/test near-geometry relationships remain:

| Retained train config | historical test | centered RMS (Å) | exact geometry | PBC match |
|---|---|---:|---|---|
| `AgTi_2p3871A_periodic_PW_PBE_force_only` | `AgTi_2p70_periodic_PW_PBE_energy_force_holdout` | 0.031609 | False | True |
| `AgTi_2p55_periodic_PW_PBE_energy_force` | `AgTi_2p70_periodic_PW_PBE_energy_force_holdout` | 0.015152 | False | True |
| `AgTi_validation_d2p60_periodic_PW_PBE_energy_force` | `AgTi_2p70_periodic_PW_PBE_energy_force_holdout` | 0.010102 | False | True |

No train/valid overlap is found under the stated screens. No historical split is silently resplit. Geometry and full-input identity are separate fields: differing PBC does not make an identical geometry independent. Persistent IDs and element identity drive comparison; missing-ID rows use an explicitly reported ordered-species fallback. Matching removes one common Cartesian translation without rotation. Exact cell/coordinate thresholds are 1e-8/1e-7 Å; near review uses all-atom RMS ≤0.05 Å/max ≤0.15 Å or separate Ag/framework RMS ≤0.15 Å. These are duplicate-review thresholds, not physical cutoffs or proof of new morphology.

## Per-frame provenance and physical limitations

| Retained train provenance | count | interpretation |
|---|---:|---|
| Verified residual PW-PBE archive | 3 | method metadata and hashes verified in this run |
| Inherited explicit method metadata | 12 | original metadata preserved; original calculation settings not independently reverified here |
| Unknown inherited method | 15 | missing per-frame method remains null/unknown; no method assigned |

All 35 train/valid/test frames have a source path/hash, source-frame index, config tag, stored method field, actual label types, atom inventory and provenance status in the manifest. Removing four unsupported LCAO records addresses one reference risk, but the **15 unknown inherited training methods still block a fully verified unified PW-PBE physical-method consistency claim**. Numerically complete labels and rank do not resolve that provenance limitation. A must review inherited reference provenance or explicitly acknowledge the uncertainty when integrating.

## Exercise and failure guards

The README contains the exact successful command, executed from `/tmp` with both explicit arguments. Independent verification reparsed all outputs, recomputed rank and hashes, checked four exclusions and exactly three residual additions, and confirmed valid/test byte preservation and blind-input isolation. Re-running on the nonempty output directory is rejected with every existing file hash unchanged. A production output directory is rejected before mutation. In-memory changed source positions are rejected, and an identical frozen input with altered PBC still registers as the same geometry. These guard checks passed; their exact results are in `exercise_checks.json`.

The generated checksum list covers the five data/manifest/diagnostic files. The publication checksum list covers the source script, README, this report, exercise checks, generated files/checksum, and all 33 source files. The checksum lists exclude themselves. If source/label validation fails, preserve the isolated run and diagnostics and mark the task failed; do not bypass validation.
