# V15 provenance-aware builder comparison

Window B job: `v15_provenance_aware_builder_revision`. Both isolated builds and independent retained-frame checks passed. Formal datasets/builders, A's tasks/logs, global manifests and training settings remain unchanged. No DFT, MACE inference/training, MD or TTM was run; V15 blind labels/summaries/outputs were not read. Only frozen blind INPUT geometries and their input manifest were used for isolation.

## Scenarios and evidence

| Check / training evidence tier | Baseline | Explicit AgTi exclusion |
|---|---:|---:|
| train / valid / test frames | 30 / 2 / 3 | 29 / 2 / 3 |
| composition rank (Ag,C,Si,Ti) | 4/4 | 4/4 |
| trace partial recoveries (`dft_method` declarations) | 14 | 14 |
| trace unresolved original method/binding | 1 | 0 |
| inherited explicit declarations, not reverified | 12 | 12 |
| newly verified residual archives | 3 | 3 |
| finite REF energies and forces | all | all |
| retained geometries / REF labels / stored metadata equal prior draft | yes | yes |
| historical valid/test bytes preserved | yes | yes |
| frozen V15 blind INPUT exact/near overlaps, all splits | 0 | 0 |
| remaining historical train–test geometry overlap pairs | 3 | 2 |

Evidence tiers apply to training rows; manifests retain individual valid/test declarations and their unverified status separately. Fourteen partial recoveries are not original calculations fully verified: original input, software/setup and convergence records remain absent. The 12 previously explicit historical rows were not upgraded to verified by this revision. Only the three residual acquisitions pass the existing archive/convergence/hash/geometry checks. Hash-bound provenance records preserve raw `source_method` and raw `dft_method` independently, declared parameter evidence (k points/cutoff/smearing where available), original-run verification status and precise evidence paths/frame/hash references.

## Label and split identity

Baseline train SHA256 is `b7edc25d671fb26e16637805b4a67493c543f66eb3302d96505dd29c34470ff7`, byte-identical to the prior isolated draft. Exclusion train SHA256 is `e358e49cf99b42d7c8ec6e482f1de5e777fa0f5da916ae9c22311ec5fc84bcf0`. The independent verifier compared every retained frame's ordered elements, positions, cell/PBC, all atom arrays including IDs/REF forces, and stored metadata/REF energy to the earlier draft. No label replacement occurs. REF/native energy and free energy remain distinct; trace alias comparisons are documented without interpreting them as independent original-run evidence.

Valid/test retain byte identity to V14 and the earlier draft. Both scenarios therefore retain historical test correlation; historical test is regression evidence, not an independent generalization score. The original V14 train–test screen had five exact/near pairs; after four LCAO exclusions it has three, and after this extra exclusion it has two. Counts are overlapping frame pairs, not unique test frames; complete pairs and thresholds are in each verification report. Frozen new V15 blind INPUT geometries have no exact/near overlap with train/valid/test under the unchanged geometry policy. This does not establish independence from parent motif lineage or unseen morphology.

## Explicit exclusion and recommendation

Both builds exclude exactly the same four force-only LCAO records and add exactly the same three verified residual labels. The optional fifth exclusion is the sole unresolved AgTi row, matched by config, source frame 13 and source file SHA256 plus trace path; exactly one match is mandatory. Config text alone cannot trigger removal. The exclusion ledger preserves exact provenance and reason; no other partially documented frame is automatically discarded.

**Recommend the 29-frame scenario for an explicitly limited numerical screening run, if A authorizes it.** It removes the sole unresolved energy/force-method binding without losing rank or any validation/test rows. Baseline remains a reproducible sensitivity comparison, not a stronger provenance choice. Exclusion reduces one risk; it does not resolve the 14 partial records, reverify the other 12 historical declarations, establish k-point convergence consistency, or make historical test independent. A owns production integration, original-record recovery/relabeling decisions, and any training authorization. No screening/training was run by B.

## Reproducibility and checks

Exact build commands are in README. `verify_scenarios.py` independently validates both runs, current source hashes, trace tier counts and exact exclusion inventory, and checks refusal of a nonempty destination and production output path. Both refusal checks exit nonzero as required. Every build refuses nonempty child directories and preserves failure diagnostics. `comparison.json` stores independent outcomes; two scenario manifests store source/output hashes, detailed ledgers, rank matrices, evidence tiers and frozen INPUT inventories. Root `SHA256SUMS.txt` covers all deliverables and source evidence without reading blind outputs.
