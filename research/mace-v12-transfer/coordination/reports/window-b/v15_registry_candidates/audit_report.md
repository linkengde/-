# v15 Ag lateral-registry candidate geometry audit (Window B)

## Scope

This publication contains label-free geometry proposals and geometric screening only. No DFT, MACE inference/training, MD, or TTM was run. The structures are new local registries derived from the existing v14 small-cluster training motifs; they do not establish independent morphology or extended-interface coverage.

The parent structures are selected by their registered `config_type` in the v14 training `train.extxyz`. Persistent `lammps_id` values, `central_pair` tags, cell, and PBC are preserved. All reference energy and force fields were removed. The full Ag sublattice is translated rigidly; all non-Ag positions remain unchanged.

## Reproducible axes and duplicate screen

The lateral axis follows the convention already used by `pbe_interface_v13_holdouts/prepare_independent_holdouts.py`: `unit(cross(r_marked_nonAg - r_marked_Ag, z_hat))`. Both signs at 0.4, 0.5, 0.6, and 0.7 Å were screened (24 translations total). This uses the repository’s established z-normal convention for its large-vacuum cluster cells.

Candidates fail the strict distance gate if any distinct-atom minimum-image distance is below 1.75 Å. Exact duplicates require matched persistent IDs, cell components within 1e-5 Å, and aligned all-atom RMS <=1e-5 Å. A near duplicate requires, after matching IDs and aligning the non-Ag framework, both Ag-registry RMS <= 0.15 Å and non-Ag framework RMS <= 0.15 Å. This checks registry similarity directly so a small Ag sublattice is not diluted by whole-structure RMS.

Every extxyz file under v13/v14 train/valid/test, v13 holdouts, v14 main holdouts, v14 parallel acquisition/holdout records, and the current v15 targeted acquisition directory was hashed and scanned. Paths, file hashes, and frame counts are in `candidate_screening_manifest.json`.

## Candidate screen

| Interface | Translation (Å) | Marked pair | Shortest total pair | Shortest Ag–nonAg pair | Source RMS / max (Å) | Nearest existing Ag-registry RMS | Decision |
|---|---:|---|---|---|---:|---:|---|
| AgC | [+0.6485, -0.2636, -0.0000] | C#604804–Ag#688433 (2.3659 Å) | C#604587–Ti#604590 (1.7413 Å) | C#604804–Ag#688433 (2.3659 Å) | 0.1429 / 0.7000 | 0.5088 Å | REJECTED BY STRICT ALL PAIR CUTOFF NOT NOMINATED |
| AgSi | [-0.6834, +0.1517, +0.0000] | Si#366934–Ag#683609 (2.5385 Å) | Ti#368832–C#368841 (1.8719 Å) | C#368645–Ag#735549 (1.9605 Å) | 0.4041 / 0.7000 | 0.7000 Å | PASS GEOMETRY SCREEN NOMINATED FOR A REVIEW |
| AgTi | [-0.0729, -0.6962, +0.0000] | Ti#479027–Ag#687247 (2.6810 Å) | Ti#480811–C#480813 (1.8722 Å) | C#480975–Ag#701905 (2.1416 Å) | 0.5292 / 0.7000 | 0.4806 Å | PASS GEOMETRY SCREEN NOMINATED FOR A REVIEW |

The candidate files are:

- `research/mace-v12-transfer/coordination/reports/window-b/v15_registry_candidates/AgC_v15_registry_candidate_label_free.extxyz` — SHA256 `ed7157c939d07a29e207f305fe371b4b7509cc771a16d519240ebbc8792d0e6d`
- `research/mace-v12-transfer/coordination/reports/window-b/v15_registry_candidates/AgSi_v15_registry_candidate_label_free.extxyz` — SHA256 `95e9980d1a2671c97865c17e3a49bda9278376d5ff550495f28802b4f84b5b70`
- `research/mace-v12-transfer/coordination/reports/window-b/v15_registry_candidates/AgTi_v15_registry_candidate_label_free.extxyz` — SHA256 `0291123e4f6db7f1f4f5e4291366048fec2fb0b7bacbf922836b16878e7b929e`

## Marked-pair nearest neighbors

### AgC

- Marked Ag atom: Ag#688433–C#604804 (2.3659 Å); Ag#688433–Ti#604600 (2.6027 Å); Ag#688433–C#604602 (2.6291 Å); Ag#688433–Ti#602709 (2.6502 Å)
- Marked framework atom: C#604804–Ti#604589 (1.9371 Å); C#604804–Ti#604796 (2.1047 Å); C#604804–Ag#688433 (2.3659 Å); C#604804–C#604602 (2.4224 Å)

### AgSi

- Marked Ag atom: Ag#683609–Si#366934 (2.5385 Å); Ag#683609–Si#367118 (2.7253 Å); Ag#683609–Ag#730721 (2.8516 Å); Ag#683609–Ag#696217 (2.8628 Å)
- Marked framework atom: Si#366934–Ag#685172 (2.0008 Å); Si#366934–Si#366936 (2.5204 Å); Si#366934–Ag#683609 (2.5385 Å); Si#366934–Ti#366935 (2.5579 Å)

### AgTi

- Marked Ag atom: Ag#687247–Ti#479027 (2.6810 Å); Ag#687247–Ag#686324 (2.8166 Å); Ag#687247–Ag#689682 (2.8194 Å); Ag#687247–Ag#685467 (2.8207 Å)
- Marked framework atom: Ti#479027–C#480805 (1.9647 Å); Ti#479027–C#479030 (2.0565 Å); Ti#479027–Ag#687247 (2.6810 Å); Ti#479027–Ti#480810 (2.7222 Å)

## Strict-cutoff finding for Ag-C

The Ag-C source already contains the non-Ag C#604587–Ti#604590 pair at 1.7413 Å. Rigid translation of Ag leaves that pair unchanged, so none of the screened Ag-C proposals can pass a literal all-atom 1.75 Å cutoff. The retained Ag-C proposal has a shortest Ag–nonAg contact above the cutoff, but is marked `REJECTED` and not nominated under the strict all-pair rule. If A intends the 1.75 Å gate to apply only to registry-involving Ag contacts, that distinction should be reviewed before any DFT assignment.

Ag-Si and Ag-Ti each have one candidate passing the strict all-pair and Ag–nonAg distance gates and the duplicate checks. They remain small-cluster registry proposals derived from existing training geometries, not independent morphologies. Window A reviews and decides any subsequent DFT labels.

## Provenance

Source dataset: `research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v14_interface_energy/data/train.extxyz`; SHA256 `3323cdbd3081fc8d91673662c5f83c8da93ba3b4d7a23428f6b70d0e7a4c9faf`.

For per-parent source-frame hashes, persistent Ag atom IDs, all 24 trial vectors and screens, nearest matched structures, complete inventory hashes, output hashes, and label-removal validation, see `candidate_screening_manifest.json`.
