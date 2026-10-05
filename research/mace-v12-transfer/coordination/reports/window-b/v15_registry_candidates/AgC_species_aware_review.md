# Ag-C species-aware framework review and registry candidate

## Scope

This is a geometry-only audit. No DFT, MACE inference/training, MD, or TTM was run. No Window A task, DFT result, or main work log was edited. B identity is `c5035b48-f43b-4da0-b8e4-e2862817f86a`.

## C#604587–Ti#604590 framework contact

C#604587 (lammps type 3) and Ti#604590 (type 1) are non-Ag atoms in the Ti₃SiC₂-derived framework. They are not the marked Ag-C registry pair (C#604804–Ag#688433). The C#604587–Ti#604590 pair recurs at short first-shell distances in eight distinct repository framework geometries sharing these persistent atom IDs:

| Representative geometry | Record class | C#604587–Ti#604590 (Å) | Median nearest Ti per framework C (Å) | C–Ti pairs <1.9 Å | Representative path/config |
|---|---|---:|---:|---:|---|
| `AgC_registry_holdout_v14_01` | DFT holdout | 1.7372 | 1.9380 | 6 | `research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v14_main_holdouts/calculations/AgC_registry_holdout_v14_01/AgC_registry_holdout_v14_01_PW_PBE.extxyz` |
| `AgC_registry_strain_acq_v14_01_periodic_PW_PBE_energy_force_v14_train` | v14 training | 1.7413 | 2.0009 | 4 | `research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v14_interface_energy/data/train.extxyz` |
| `AgC_secondary_2.0A_LCAO` | v11 training | 1.7660 | 2.0025 | 6 | `research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v11_energyfocus/data/train.extxyz` |
| `AgC_secondary_2.4A_LCAO` | v11 training | 1.7660 | 1.9821 | 7 | `research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v11_energyfocus/data/train.extxyz` |
| `AgC_validation_d2p30_periodic_PW_PBE_energy_force` | v13 validation data | 1.7660 | 1.9821 | 6 | `research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v13_interface_energy/data/train.extxyz` |
| `AgC_lateral_registry_holdout_v13_periodic_PW_PBE_energy_force_v14_train` | v14 training | 1.7823 | 1.9767 | 6 | `research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v14_interface_energy/data/train.extxyz` |
| `AgC_cluster_PW_PBE` | DFT cluster output | 1.7660 | 1.9821 | 6 | `research/mace-v12-transfer/periodic_interface_v4/pbe_cluster_transfer_pw_20261002_run3_mixer_retry2/AgC_cluster_PW_PBE.extxyz` |
| `AgC_residual_shell_v15_01` | DFT residual target | 1.8102 | 1.9981 | 5 | `research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v15_targeted_acquisition/calculations/AgC_residual_shell_v15_01/AgC_residual_shell_v15_01_PW_PBE.extxyz` |

In the candidate, C#604587’s nearest Ti shell is Ti#604590 1.7413 Å, Ti#604584 1.8395 Å, Ti#604586 1.9322 Å, Ti#604589 2.0886 Å. Ti#604590’s nearest C shell is C#604587 1.7413 Å, C#604596 1.8963 Å, C#604602 2.0229 Å, C#602563 2.0277 Å. This is a recurring framework first-shell adjacency, not an isolated Ag-induced contact. Its distance is at the short/compressed end of this motif’s repeated DFT/data examples (1.7372–1.8102 Å), while per-carbon nearest-Ti medians are about 1.94–2.00 Å. The best-supported description is a recurrent, locally compressed/distorted framework C–Ti bond; the repository does not establish it as an ideal bulk equilibrium bond. A uniform 1.75 Å all-pair cutoff therefore cannot chemically reject this framework pair.

## Existing Ag-C geometry search

Scanned 52 existing `.extxyz`/`.xyz` files (163 frames) under `research/mace-v12-transfer`, excluding only the new output candidate itself. This inventory includes the earlier B label-free proposal. It contains 36 stored Ag-C central-pair frame records representing 13 exact coordinate geometries. The old B proposal is one of those geometries; its coordinates are not new in this task. The search included v11–v14 train/validation/test records, calculation inputs and DFT outputs, v13/v14 holdouts/acquisitions, and B’s existing proposals. Of 28 records matching the candidate’s formula and all persistent IDs, 27 are existing dataset or calculated-DFT geometry references after excluding B’s prior proposal. These other geometries are represented in train/split records or have calculated DFT records.

### Candidate

The structure is published as `AgC_v15_registry_candidate_species_aware.extxyz`. It is the prior B 0.7000 Å rigid translation, now reassessed with an element-aware framework check; its coordinates are exactly unchanged from `AgC_v15_registry_candidate_label_free.extxyz` (maximum MIC coordinate difference 0.0 Å). It derives from the v14 training frame `AgC_registry_strain_acq_v14_01_periodic_PW_PBE_energy_force_v14_train` in `research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v14_interface_energy/data/train.extxyz` (source file SHA256 `3323cdbd3081fc8d91673662c5f83c8da93ba3b4d7a23428f6b70d0e7a4c9faf`; source-frame geometry SHA256 `72d649542158aade39a3cf50dc73759ad25df8a7ff7b2b117f1aa0c6ffc513a3`). Ag IDs 688433 and 751743 are translated by `[+0.64849057, -0.26355263, 0.00000000] Å` (magnitude 0.7000 Å); non-Ag framework coordinates are unchanged from the source. The structure has 48 atoms (C13Ag2Si10Ti23), preserves all persistent IDs and the marked pair C#604804–Ag#688433, cell lengths 28.1400 Å, and PBC `[True, True, True]`. It contains no energy or force labels.

The species-aware registry screen compares each Ag-containing pair type with the shortest distance observed in existing labeled DFT geometries of this same small-cluster motif. These empirical values are reference observations, not universal chemical cutoffs:

| Ag pair type | Candidate minimum | Shortest same-motif labeled DFT minimum | Reference config | Result |
|---|---:|---:|---|---|
| Ag-C | 2.3659 Å | 2.0000 Å | `AgC_d2p0_periodic_PW_PBE_energy_force` | PASS: not shorter than reference |
| Ag-Ti | 2.6027 Å | 2.3069 Å | `AgC_residual_shell_v15_01` | PASS: not shorter than reference |
| Ag-Si | 3.9473 Å | 3.7884 Å | `AgC_registry_holdout_v14_01` | PASS: not shorter than reference |
| Ag-Ag | 2.9182 Å | 2.8457 Å | `AgC_d2p4_periodic_PW_PBE_energy_force` | PASS: not shorter than reference |

Duplicate comparison matched persistent IDs and species against 27 same-formula Ag-C data/DFT frame records. Exact geometry hits: 0. Registry-near hits: 0, using Ag-only RMS ≤0.15 Å and non-Ag framework RMS ≤0.15 Å after common Cartesian translation alignment of these centered large-vacuum clusters (exact duplicate tolerance 1e-5 Å, plus matching cell and PBC). The nearest reference is `AgC_validation_d2p30_periodic_PW_PBE_energy_force` in `research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v13_interface_energy/data/train.extxyz` frame 23: Ag-registry RMS 0.5088 Å, framework RMS 0.0352 Å, all-atom RMS 0.1094 Å. Whole-structure RMS is diluted when only two Ag move, so the Ag-registry RMS is the relevant distinctness measure.

**Disposition:** passes the empirical species-pair contact screen and is not duplicated by the existing training/validation/test or calculated DFT configurations searched. It is a valid new local Ag registry relative to labeled records, but it reuses a small-cluster framework and is the same geometry as B’s earlier unlabeled proposal. It is not an independent morphology or extended-interface holdout. Window A decides whether to freeze it and assign DFT.

## Provenance and hashes

Candidate SHA256: `f528a051246a1a89546da6f55ee3e3b66222cd08caefca039d77989049d25135`. Input file hashes, all Ag-C central-pair records, the 27 duplicate comparisons, framework C–Ti references, pairwise screen, and label-exclusion checks are in `AgC_species_aware_manifest.json`. The flat checksum list is `AgC_species_aware_SHA256SUMS.txt`.
