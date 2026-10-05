# Window B v15 dataset integrity and label provenance audit

## Scope and recommendation

Read-only review of committed v11–v14 datasets/builders, v15 residual acquisition records, and the three frozen v15 registry inputs. No DFT, MACE, MD/TTM, dataset generation, or builder execution occurred. No A temporary checkpoint or live process was inspected. Audit inputs were read from `main` commit `56e4657e3bf9c0ea903c87963658c3a433dfe74a`; their SHA256 values are in the companion JSON and checksum list.

**Do not mark v15 ready to train yet.** All three residual labels have compact verified archives on `main`; A's task-state values for the three fresh registry labels are recorded in JSON. Holdout label contents were not opened. The v15 builder has not been produced. Before the PW-PBE objective is built, exclude or replace four legacy force-only LCAO records and preserve all three fresh registry labels as blind post-freeze tests.

## Split lineage, provenance and label counts

| Dataset | Train frames / energy / force | Valid frames / energy / force | Test frames / energy / force | Result |
|---|---:|---:|---:|---|
| v12 | 23 / 19 / 23 | 2 / 2 / 2 | 3 / 3 / 3 | counts and output hashes match manifest |
| v13 | 26 / 22 / 26 | 2 / 2 / 2 | 3 / 3 / 3 | counts and output hashes match manifest |
| v14 | 31 / 27 / 31 | 2 / 2 / 2 | 3 / 3 / 3 | counts and output hashes match manifest |

The JSON companion lists **every frame** in v12, v13 and v14 train/valid/test with config tag, source-method field as stored, label keys, formula, ID count, cell/PBC, geometry hash, source path and file SHA256. In v14 train, 12 energy+force frames explicitly record the PW-PBE method; 15 inherited energy+force frames lack per-frame `source_method`, and all 4 force-only frames also lack it. Newer builder-added frames record `GPAW 26.7.0 PW-PBE 500 eV Gamma single point`; do not silently infer inherited-frame method from a later frame.

The four remaining force-only LCAO training frames are:

| v14 train frame | atoms/formula | reference labels in extxyz | per-frame method field | PBC |
|---|---|---|---|---|
| `AgC_primary_2.0A_LCAO` (frame 9) | 48 / `C12Ag13SiTi22` | energy: no; forces: yes | missing; LCAO appears only in config tag | `[False, False, False]` |
| `AgC_primary_2.4A_LCAO` (frame 10) | 48 / `C12Ag13SiTi22` | energy: no; forces: yes | missing; LCAO appears only in config tag | `[False, False, False]` |
| `AgTi_2p1A_LCAO` (frame 11) | 49 / `C11Ag28SiTi9` | energy: no; forces: yes | missing; LCAO appears only in config tag | `[False, False, False]` |
| `AgTi_2p7A_LCAO` (frame 12) | 49 / `C11Ag28SiTi9` | energy: no; forces: yes | missing; LCAO appears only in config tag | `[False, False, False]` |

The v12 builder replaces `AgSi_2p2A_LCAO`, `AgSi_2p8A_LCAO`, and `AgC_secondary_2.0A_LCAO`/`AgC_secondary_2.4A_LCAO`; those are different frames from the four remaining `AgC_primary_*`/`AgTi_*` records above. The v11 manifest's global note says “All are PBE/500 eV GPAW labels,” but it does not supply per-frame source logs or resolve the explicit LCAO tags and missing-energy fields. The committed repository has no per-frame original LCAO inputs, basis settings, GPAW version, XC/cutoff/k-point/smearing settings, convergence record, or original output hashes for these four forces. Do not infer those settings.

The v12/v13 frame audits show the trainer's unmasked energy metric included the four missing energy values as zeros (1800.8 and 1697.4 meV/atom, respectively); only 19/22 frames had real energy labels. A future v15 objective must remove these rows or mask missing energy targets; zero-filled values are not labels.

## Leakage and independence screen

Comparison maps persistent `lammps_id` values and requires the same element per ID. For centered large-vacuum clusters, one common Cartesian translation is removed. Exact duplicate means same cell within 1e-8 Å, same PBC and same-ID RMS/max ≤1e-8 Å. Near screen means same cell and same-ID RMS ≤0.05 Å and maximum ≤0.15 Å; PBC is reported separately. For the new registry comparisons, also require both Ag-sublattice RMS and framework RMS ≤0.15 Å after framework alignment, so the moving Ag atoms are not diluted by whole-system RMS.

The V14 frozen AgTi test `AgTi_2p70_periodic_PW_PBE_energy_force_holdout` has the same 49 IDs, elements, cell and Cartesian coordinates as train frame `AgTi_2p7A_LCAO` (RMS 0, max 0), but PBC differs: train `[False,False,False]`, test `[True,True,True]`. It is therefore a train/test **coordinate leak** even though the complete calculator inputs differ. The train/test force-vector disagreement at those coordinates is 0.0852 eV/Å RMS, max atom-vector difference 0.2013 eV/Å. The test is also near the PW-PBE train points `AgTi_2p55_periodic_PW_PBE_energy_force` (RMS 0.0152 Å, max 0.0750 Å) and `AgTi_validation_d2p60_periodic_PW_PBE_energy_force` (RMS 0.0101 Å, max 0.0500 Å). Do not count this frozen test as independent geometry evidence. These split structures are unchanged in v12–v14, so repeated version scores are not independent samples.

The prior v14 **registry holdout** `AgTi_registry_holdout_v14_01` is a different record from the leaked 2.70 Å frozen test; the force-separation plan reports its nearest training coordinate RMS 0.4153 Å/max 0.5746 Å. It was already scored and used for diagnosis, so retain it as historical regression evidence, not as a fresh v15 blind test.

The three new v15 registry inputs have no labels, match their SHA manifest and have zero exact/near duplicate matches to v14 training under the registry threshold. The 0.05 Å RMS and 0.15 Å maximum displacement values above are a review threshold for duplicate triage, not a universal chemistry cutoff. Nearest same-ID v14 train references are:

| Frozen input | Nearest training config | Ag registry RMS (Å) | framework RMS (Å) | A status (label record; current job; published archive) |
|---|---|---:|---:|---|
| `AgC_registry_holdout_v15_01` | `AgC_registry_strain_acq_v14_01_periodic_PW_PBE_energy_force_v14_train` | 0.7000 | 0.0000 | `queued`; job `running; SCF iteration 0; elapsed 48s; updated 2026-10-05T06:00:43.778533+00:00`; archive paths `0`; contents not inspected |
| `AgSi_registry_holdout_v15_01` | `AgSi_registry_strain_acq_v14_01_periodic_PW_PBE_energy_force_v14_train` | 0.7000 | 0.0000 | `queued`; job `queued; SCF iteration 0; updated 2026-10-05T05:59:55.426433+00:00`; archive paths `0`; contents not inspected |
| `AgTi_registry_holdout_v15_01` | `AgTi_validation_d2p60_periodic_PW_PBE_energy_force` | 0.4810 | 0.0176 | `queued`; job `queued; SCF iteration 0; updated 2026-10-05T05:59:56.965331+00:00`; archive paths `0`; contents not inspected |

By contrast, the three residual acquisition inputs are intentionally close to their v14 training parents (the JSON gives per-interface Ag/framework RMS). They are training candidates, not independent validation. At the analyzed main commit, all three residuals have compact verification `PASS`, converged four-rank summaries and matching finite label arrays. Holdout task states are shown above/JSON; their label contents were deliberately not inspected and must stay outside v15 train, validation, checkpoint selection, and hyperparameter selection until the model is frozen.

## Expected v15 counts and composition support

The v14 train split has 31 structures, 27 energy labels and 31 force arrays; four structures are the force-only LCAO records above. If all three verified PW-PBE residual labels enter training **and the four unsupported LCAO frames are excluded**, the expected v15 training set is 30 structures with 30 energy and 30 force labels. The currently archived residual labels support that projected count. If the LCAO rows were retained, the naive count would be 34 structures, 30 energy labels and 34 force arrays, mixing four unverified force references into the PW-PBE objective. The three frozen registry holdouts are separate, label-pending evaluations and are not included in those counts.

The v14 27-frame energy composition matrix has rank 4 in columns `[Ag,C,Si,Ti]`. The V14 training command used `energy_weight=100` and `forces_weight=1000` (`run_v14_training.sh`); hold these fixed for the first controlled v15 comparison. Adding the three residual compositions projects to 30 energy rows and rank 4. This is a composition-identifiability check only; it does not validate the labels or model. The exact rank and singular values are in JSON.

## Builder acceptance checklist

1. Add only verified compact PW-PBE archives with input/output/summary hashes matching manifests, convergence and four MPI ranks; require finite energy and an `(N,3)` finite force array.
2. Exclude the four legacy force-only LCAO frames unless matched PW-PBE labels and complete provenance are supplied. Never encode missing energies as zero targets.
3. Emit frame-level provenance: config type, source path/SHA, source method, label fields, element-by-ID mapping, cell and PBC. Reject unexplained ID/order or cell changes.
4. Check exact and near duplicates by ID-matched species/coordinates, cell and PBC across every split. Flag AgTi 2.70 Å as historical leakage; do not use it as independent evidence.
5. Keep all three frozen registry inputs/labels outside training, validation, checkpoint and hyperparameter choices. Score only after checkpoint freeze; keep previous v14-scored structures as historical regressions.
6. Assert exact frame/energy/force counts and full-rank energy composition (4/4); distinguish energy-labeled rows from force-only rows and report singular values.
7. Preserve v14 frozen split bytes when carried forward, and verify/hash every v15 input and generated split against its manifest.
8. Keep the first v15 training comparison at the V14 baseline settings; use the new holdouts only for post-freeze scoring.

## Audit disposition

| Item | Verdict |
|---|---|
| v12-v14 split counts, manifest hashes and byte lineage | PASS |
| Four legacy LCAO force-only references for a unified PW-PBE loss | FAIL pending matched labels/provenance; exclude for now |
| V14 AgTi frozen test as independent geometry evidence | FAIL; same coordinate set as LCAO train plus two near PW-PBE train points |
| AgC/AgSi/AgTi residual archives | PASS for archive/label integrity after A's builder checks |
| New v15 holdout geometry uniqueness and frozen hashes | PASS at input stage; label content kept uninspected |
| Projected composition rank | PASS (4/4) |
| v15 training readiness | UNKNOWN / NOT READY; no v15 builder output; holdouts remain reserved and uninspected |

Input file hashes, all 95 v12-v14 split-frame records, duplicate findings, label archive evidence, and projected counts are in `v15_dataset_integrity_audit.json`. At this commit, the A task's `labels` map lists the three v15 holdouts as queued, while `current_jobs` reports AgC as running; no holdout archive paths are committed. This bookkeeping discrepancy is preserved in JSON. No holdout label files were opened. The checksum file covers the audited inputs plus report and JSON; it does not self-hash.
