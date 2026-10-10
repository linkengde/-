# Ti₃SiC₂ backbone force-response DFT batch — B assignment

## Purpose

Label the four source-audited, fixed-cell Ti₃SiC₂ backbone response geometries prepared in the preceding B report. The selected V18 model has no pure Ti₃SiC₂ training frame and gives a large force error on the public bulk control. These DFT results will provide candidate backbone-response data; they are not independent validation, and they must not be integrated into training until A has defined the data split and verified the complete model-entry workflow.

## Inputs and exact hashes

Source report: `research/mace-v12-transfer/coordination/reports/window-b/v19_Ti3SiC2_backbone_and_interface_acquisition_design_window_b/`.
Its SHA-256 inventory was verified before this task was assigned. Four inputs are listed in `proposal_manifest.json`; their content hashes are:

| Geometry | New calculation label | Input SHA-256 |
|---|---|---|
| `Ti3SiC2_Ti4f_id0003_z_m020A_PROPOSAL.extxyz` | `Ti3SiC2_Ti4f_id0003_z_m020A_k10x10x4_v19` | `5529a040f0f35957f826b78754910970741bff741474e321dc1646a5c0c12bbf` |
| `Ti3SiC2_Ti4f_id0003_z_p020A_PROPOSAL.extxyz` | `Ti3SiC2_Ti4f_id0003_z_p020A_k10x10x4_v19` | `5829ebcbb9ae331f6cf1dd7dc5abca0787c07667a55778e3c1eea241938b227b` |
| `Ti3SiC2_C4f_id0009_z_m020A_PROPOSAL.extxyz` | `Ti3SiC2_C4f_id0009_z_m020A_k10x10x4_v19` | `19f1a4011620f5c7f245d0a8c893964aa8ebd1eb647f8a9b6f4547a75cac586f` |
| `Ti3SiC2_C4f_id0009_z_p020A_PROPOSAL.extxyz` | `Ti3SiC2_C4f_id0009_z_p020A_k10x10x4_v19` | `6308961f33ad5fae3f58930570a47044a236e062de243821c20577c78955ff72` |

The geometry generator preserved atom IDs/order, cell and PBC. Each proposal moves exactly one representative 4f atom by ±0.02 Å along z from the verified 48-atom Ti₃SiC₂ parent. The preparation report's overlap, uniqueness, identity and source-package hash checks passed. These four input hashes are not present in the existing pure-phase calculation manifest or archive labels.

## Registered method

Use the exact verified method of `Ti3SiC2_baseline_k10x10x4_v19`, including PBE, 500 eV plane-wave cutoff, 0.10 eV Fermi–Dirac smearing, 10×10×4 k-points, three-dimensional periodicity, 0.1 eV smearing, the existing mixer and `1e-5` SCF thresholds, GPAW 26.7.0, ASE 3.29.0, and gpaw-data 1.2.1. The method manifest hash is `702df71e586a61cd5785afc62653d7135a707032c91b72842b035b62198fd709`. Keep native extrapolated energy and force-consistent free energy separate.

## B execution sequence

1. Sync `main`; independently verify the source report inventory and each input SHA above. Confirm the baseline archive remains verified and none of the four new labels, output directories, or launcher logs already exists.
2. Register only these four exact inputs in `pbe_interface_v19_pure_phase_controls`, using the baseline method unchanged, `owner_task=window-b`, B's registered instance ID, and the four exact labels. Add the labels to B's task card and update the package inventory. Publish this registration before launching any job.
3. Run the four exact labels sequentially, in the table order, through the existing `run_queue.py`. Use its four-rank/one-thread preflight, progress publication, archive verifier and publisher. Do not overlap runs. Stop on the first failed or unconverged job and preserve its checkpoint; do not restart a label or reuse an existing output.
4. Keep `.gpw` checkpoints local and out of Git. After each verified archive, check disk reserve before proceeding. A completed checkpoint may be removed only when its GPAW process has exited, its compact archive passes verification, and disk reserve requires cleanup; never remove an active or failed-run checkpoint.
5. Publish a report with convergence status, archive verification and hashes for every label. For each displaced atom, report the ± energy changes relative to the unchanged parent and the symmetric z-force response `-(Fz(+0.02)-Fz(-0.02))/0.04` in eV/Å², plus the all-atom force-vector RMS/max against the parent. Keep native/free-energy conventions separate.

## Scope limits

This batch is fixed-cell bulk backbone response only. Do not relax the cell or atoms, do not run interface points (the COD source bytes remain unverified and the slab numerical recipe is pending A's review), and do not train/infer a model, access dev/test/holdout/sealed labels, integrate the results into a dataset, edit A-owned files, or rerun the completed baseline. The four labels share one parent and site family; they cannot serve as independent model validation or establish overall Ti₃SiC₂ force accuracy.
