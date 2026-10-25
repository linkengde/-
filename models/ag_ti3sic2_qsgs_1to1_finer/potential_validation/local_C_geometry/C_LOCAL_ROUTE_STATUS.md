# C route: local geometry prepared; interaction-model gate blocked

Date: 2026-10-25. This directory contains a source-centered local crop from the existing QSGS model and an audit of repository-supported Ag–Ti/Si/C interaction options. The crop is a geometry artifact only; it is not a relaxed state and has not been used for dynamics.

## Local geometry

`build_local_C_crop.py` extracts a 160 × 160 Å lateral region centered beneath the existing Gaussian heat-source center, retains the complete 159.392 Å material thickness and 100 Å top vacuum, and uses x/y nonperiodic cut boundaries. It preserves source atom positions and records a source-ID map. The planned future local-cell boundary has a 20 Å fixed side buffer and 20 Å fixed bottom layer; those constraints are not yet dynamically validated.

The crop contains 252,123 atoms: 59,046 Ag and 193,077 Ti/Si/C atoms, with local Ag:ceramic mass ratio 1.0145:1. The nominal source core contains 1,275 Ag atoms; a wider r ≤ 60 Å, depth ≤ 100 Å buffer contains 15,080 Ag atoms. The nearest Ag-to-skeleton distance is 2.333 Å; no pair is closer than 2.2 Å in the extracted coordinates. Within the source core, the nearest distances are 2.334 Å to Ti, 2.338 Å to Si, and 2.342 Å to C. These are geometric distances from the as-built crop, not relaxed bond lengths. The largest local connectivity components contain 58,764/59,046 Ag atoms (99.52%) and 192,870/193,077 skeleton atoms (99.89%); the remaining small components are retained and reported, not silently deleted.

The LAMMPS geometry-parse check reads the compressed data file and confirms 96,028 Ti, 32,087 Si, 64,962 C, 59,046 Ag, and 126,893 atoms in the planned inner block. It uses `pair_style zero`, `run 0`, and no dynamics; it verifies file syntax and counts only.

## Ag–skeleton interaction audit

No interaction option presently in the repository qualifies for the planned cold-state or TTM-MD gate:

| Candidate | Repository evidence | Decision for this local model |
|---|---|---|
| Morse calibration-v0 | The include itself says “NOT YET A PRODUCTION POTENTIAL”; its starting values come from a single-plane HCP surrogate of reported adhesion targets. Existing QSGS tests show large interface-force contributions and failed cold baselines. | Do not use as a validated cross potential. |
| DFT interface labels | Existing PW-PBE records provide selected fixed geometries. The separate Ag–C, Ag–Si, and Ag–Ti distance-validation inputs each displace one central pair to a single target distance (2.30, 2.60, and 2.60 Å respectively); these do not form matched separation curves sufficient to uniquely fit and validate a transferable Ag–X pair function. | Useful source data, insufficient by themselves for a production cross potential. |
| MACE v14–v19 | v14, v15, and v16 fail force/separation-force screening; v17 and v18 fail development evaluation; v19 explicitly fails force-vector and separation-force development gates and is not qualified for production or impact simulation. The v19 development set is three structures from one source family, and slab numerical convergence remains unresolved. | Do not use for MD/TTM or claim cold-state acceptance. |
| WCA-only B obstacle | It prevents selected overlaps but has no Ag–skeleton attraction. Whole-layer, ROI, and source-core cold checks all continue to show structural adjustment. | Retain only as a failed exploratory geometry control. |

The local crop is therefore ready for potential evaluation but not for an interaction-bearing `run 0`, minimization, or 300 K dynamics. Running a zero-force or WCA-only test on it would not answer the C-route physical question.

Repository evidence used in this audit: [`potential.interface_calibration_v0.inc`](../../potential.interface_calibration_v0.inc); [`v14 validation assessment`](../../../../research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v14_interface_energy/results/v14_validation_assessment.json); [`v15 evaluation`](../../../../research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v15_interface_energy/evaluation_01/evaluation.json); [`v16 evaluation`](../../../../research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v16_interface_energy/evaluation_01/evaluation.json); [`v17 development evaluation`](../../../../research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v17_interface_energy/evaluation_dev_01/evaluation.json); [`v18 development evaluation`](../../../../research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/evaluation_dev_01/evaluation.json); [`v19 candidate status`](../../../../research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v19_residual_response/README.md); and [`DFT distance-scan manifest`](../../../../research/mace-v12-transfer/periodic_interface_v4/pbe_interface_energy_additions_v12/validation_distance_input_manifest.json).

## Gate and next concrete work

Before local relaxation, obtain one of the following: (1) a force-qualified Ag–Ti/Si/C model with evidence in matching local environments, or (2) a versioned effective interface model fitted to a matched set of DFT distance/registry/force data and accompanied by a narrow sensitivity range. The existing records do not currently supply either. Once selected, apply that model to this crop, perform one finite local reconstruction, then a 300 K NVT-to-NVE cold gate. No extra FIRE/NVT cycles or heating should be applied to this coordinate crop until the interaction-model gate is met.

## Reproduction

From the finer-model directory:

```bash
python potential_validation/local_C_geometry/build_local_C_crop.py
OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp \
  -in potential_validation/local_C_geometry/in.geometry_parse_check \
  -log potential_validation/local_C_geometry/log.geometry_parse_check
```

The geometry-parse input intentionally uses zero interactions and `run 0`; it does not test forces, stability, melting, or TTM.
