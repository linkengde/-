# v15 blind-holdout geometry design audit (Window B)

## Scope and method

This is a read-only geometry audit of the records present at `e465269` plus the B progress publication. It used the v14 per-atom residual report, the v14 validation assessment, the v15 acquisition input manifest and coordinates, and the v12/v13/v14 train/valid/test and v13/v14 holdout/acquisition geometries. No DFT, MACE inference/training, MD/TTM, geometry generation, or edits to dataset/holdout inputs were performed.

Candidate-to-training comparisons match atoms by persistent `lammps_id` and element in the same periodic cell, then measure minimum-image Cartesian displacement. The v14 holdouts use slightly different cubic vacuum cells from the v15 acquisition inputs; for those comparisons the audit records cell lengths, persistent pair IDs, contact distances, and neighbor environments instead of treating raw whole-cell RMS as an independence test.

The v14 model evaluated in the residual report has SHA256 `41bd411355b034b0dab85dbdbda932ecfc394122cb67a2b51a550bd04207645d`. The screen remains `SCREEN_FAIL`.

## v15 acquisition lineage and overlap with v14 training

All three proposed v15 acquisition inputs are small local perturbations of configurations that already appear in the v14 training split. The rows below are same-ID, same-cell comparisons to the corresponding labeled v14 training frame; their small nonzero displacement comes from the v15 local perturbation.

| Interface | v15 acquisition input and source | Corresponding v14 training configuration | Atom displacement from that training frame (RMS / max) | Marked pair and input distance |
|---|---|---|---:|---|
| Ag-C | `AgC_residual_shell_v15_01`, from `AgC_registry_strain_acq_v14_01` | `AgC_registry_strain_acq_v14_01_periodic_PW_PBE_energy_force_v14_train` | 0.0288 / 0.0602 Å | C#604804–Ag#688433, 2.2571 Å |
| Ag-Si | `AgSi_residual_shell_v15_01`, from `AgSi_registry_strain_acq_v14_01` | `AgSi_registry_strain_acq_v14_01_periodic_PW_PBE_energy_force_v14_train` | 0.0520 / 0.1270 Å | Si#366934–Ag#683609, 2.4505 Å |
| Ag-Ti | `AgTi_residual_shell_v15_01`, from `AgTi_registry_probe_v13_01` | `AgTi_registry_probe_v13_01_periodic_PW_PBE_energy_force_v14_train` | 0.0440 / 0.0886 Å | Ti#479027–Ag#687247, 2.6010 Å |

These are useful acquisition perturbations, but they are correlated with already-trained geometries and must not be counted as v15 blind tests.

The marked-pair first shells in the v15 inputs are:

- **Ag-C:** C#604804 has Ti#604589 at 1.934 Å and Ti#604796 at 2.102 Å; its Ag#688433 contact is 2.257 Å. Ag#688433 also neighbors C#604602 at 2.244 Å and Ti#604600 at 2.307 Å.
- **Ag-Si:** Si#366934 contacts Ag#683609 at 2.451 Å and neighbors Ti#366935 at 2.472 Å, Si#366936 at 2.568 Å, and Ag#685172 at 2.579 Å. Ag#683609 neighbors Ti#367117 at 2.540 Å and Si#366936 at 2.746 Å.
- **Ag-Ti:** Ti#479027 contacts Ag#687247 at 2.601 Å and neighbors C#480805 at 2.001 Å, C#479030 at 2.117 Å, and Ag#682693 at 2.325 Å. Ag#687247's next Ag neighbors are #689682 at 2.713 Å, #686324 at 2.736 Å, and #716228 at 2.803 Å.

## Coverage of the localized v14 residuals

The v15 generator perturbs atoms within 4.5 Å of its listed focus IDs. The following is geometric coverage only; it does not predict the DFT forces of the pending acquisition labels.

| Interface / v14 holdout | Largest localized v14 errors | Distance from those atoms to the v15 focus set | Coverage assessment |
|---|---|---|---|
| Ag-C / `AgC_registry_holdout_v14_01` | Ti#604590: 0.3241 eV/Å; Ti#604806: 0.2978; C#604602: 0.2875 | 3.260, 3.690, 2.133 Å; none is an explicit focus ID | All three lie inside the perturbed shell, so the acquisition samples the same local region. It remains a perturbation of the v14 training motif rather than an independent holdout. |
| Ag-Si / `AgSi_registry_holdout_v14_01` | Si#366934: 0.2793 eV/Å; C#366931: 0.2359; Ti#368832: 0.2249 | 0.000 Å (explicit focus), 2.657 Å, 5.217 Å | Partial. The Si and C sites are in the perturbed shell; Ti#368832 is outside the 4.5 Å radius and its local environment is not directly perturbed by this input. |
| Ag-Ti / `AgTi_registry_holdout_v14_01` | C#480805: 0.4652 eV/Å; Ti#479027: 0.3022; Ti#480812: 0.3002 | 2.001 Å, 0.000 Å (explicit focus), 0.000 Å (explicit focus) | All three lie in the perturbed shell; this is good local coverage of the reported Ag-Ti residual environment, but not independent motif coverage. |

The corresponding v14 holdout marked-pair distances are 2.2491 Å (Ag-C), 2.4387 Å (Ag-Si), and 2.5637 Å (Ag-Ti). The holdouts share the same persistent pair IDs and first-shell species as the acquisition inputs. Their recorded families are `AgC_existing_cluster_motif`, `AgSi_existing_cluster_motif`, and `AgTi_existing_cluster_motif`; their cubic vacuum cells differ slightly (27.72, 28.28, and 27.72 Å, respectively). The pairs and environments are therefore related small-cluster motifs, not evidence of independently sourced morphologies.

## Fresh blind-holdout availability

No defensible, not-yet-used distinct registry or motif for a fresh blind geometry per interface is present in the inspected v4 records:

- **Ag-C:** the v15 source is the v14 strain-acquisition geometry already in v14 training. The older `AgC_cluster_PW_PBE.extxyz` source aligns by persistent IDs to v14 training frame `AgC_validation_d2p30_periodic_PW_PBE_energy_force` (RMS 0.0099 Å, max 0.0487 Å). The separate `AgC_primary_2.0A_LCAO`/`AgC_primary_2.4A_LCAO` geometries are also in v14 training. The v13 lateral-registry holdout was added to v14 training, and the v14 Ag-C holdout has already been scored and used in the residual diagnosis.
- **Ag-Si:** the v15 source is the v14 strain-acquisition geometry already in v14 training. The older v8 `AgSi_cluster_PW_PBE.extxyz` source aligns to v14 training frame `AgSi_validation_d2p60_periodic_PW_PBE_energy_force` (RMS 0.0144 Å, max 0.0729 Å). The v13 lateral-registry holdout was added to v14 training; the v14 Ag-Si holdout has already been scored and localized.
- **Ag-Ti:** the v15 probe source is in v14 training. The d2.55 and `AgTi_validation_d2p60` geometries are in the v14 training split, and the d2.70 geometry is in the frozen v14 test split. The v14 Ag-Ti holdout has already been scored and used for residual localization. No separate unscored Ag-Ti registry/motif input appears in the inspected records.

Partition check: v14 contains 31 training, 2 validation, and 3 frozen-test frames. The exact v15 source frames and the Ag-C/Ag-Si v13 lateral-registry records are in training; the available C/Si gap tests and AgTi d2.70 test have already been scored. The older Ag-Si distance scans are likewise in training. No unused geometry with a distinct registry and provenance remains among the inspected records.

The other gap-series structures (`C_gap`, `Si_gap`, and `Ti_gap`) are already assigned to v14 train/validation/test and do not supply a fresh blind example for these marked Ag-X contact motifs. Reusing a previously scored v14 holdout could keep it out of v15 training, but it would not be a fresh blind geometry for this error-guided cycle.

**Recommendation:** report a source-geometry gap for all three interfaces; do not synthesize or nominate coordinates from the existing files as fresh blind structures. Window A should supply one new geometry per interface from a distinct registry or independently sourced motif. Freeze and hash those inputs before any v15 scoring, and keep them out of v15 training and validation until their blind scores are recorded. If new source geometries cannot be obtained, document that limitation rather than claiming fresh independent coverage.

## Reproducible source paths

- `coordination/reports/window-b/v14_force_localization.json`
- `periodic_interface_v4/mace_periodic_v14_interface_energy/results/v14_validation_assessment.json`
- `periodic_interface_v4/pbe_interface_v15_targeted_acquisition/input_manifest.json`
- `periodic_interface_v4/mace_periodic_v14_interface_energy/data/{train,valid,test}.extxyz`
- `periodic_interface_v4/pbe_interface_v13_holdouts/inputs/`
- `periodic_interface_v4/pbe_interface_v14_main_holdouts/inputs/`
- `periodic_interface_v4/pbe_interface_v14_parallel_acquisition/inputs/`
