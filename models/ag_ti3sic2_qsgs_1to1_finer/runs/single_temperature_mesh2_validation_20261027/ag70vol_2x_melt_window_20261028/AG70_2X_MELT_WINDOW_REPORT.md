# Ag 70:30 vol% 2× melt-window validation

## Status

The 2× restart checkpoint at 0.35 ps passed integrity and physics-identity validation. A same-checkpoint 0.01 ps baseline and a profiling run completed, and an implementation-only optimization reproduced all saved state arrays bit-for-bit. The long continuation is running from the validated 0.36 ps checkpoint in heating_0p36_to_pool/. Its first complete 0.05 ps segment reached 0.41 ps, passed restore-only validation, and remained fully solid (peak 790.581 K); no 2× melting result is claimed yet.

## Fixed model conditions

- Ag/Ti₃SiC₂ composition: 70:30 by volume, using the existing QSGS geometry.
- Refinement: 2×, mesh 200 × 200 × 100.
- Case A, existing source, thermophysical candidates, phase-change rules and time-step controller are unchanged.
- External faces remain insulated; Ti₃SiC₂ has zero storage and conductivity in Case A.
- The 2300 K protection and 0.30 ps post-pool follow-up are enabled. The configured time cap is 5.0 ps; the driver stops earlier after its pool follow-up condition or the temperature guard.
- No 1× simulation or QSGS generation was repeated.

## Checkpoint and benchmark verification

The supplied checkpoint ag70vol_resume_2x/checkpoints/caseA_ag70_2x_heat_0p35ps.npz restored at 0.35 ps with peak temperature 777.1887768 K, zero Ag liquid fraction, injected energy 5.429756311829655e-16 J and energy residual 2.3271396704019848e-29 J (relative residual 4.29e-14). All four 200 × 200 × 100 fields were finite, and the embedded geometry identity matched SHA-256 03ef47b1dd534b3215cf03f42c8cf7e887dcdfcc185b33bab7453cd838c1faac.

| Run | Wall time | Peak RSS | End state | Energy residual |
| --- | ---: | ---: | --- | ---: |
| Baseline, 0.35→0.36 ps | 287.850 s | 781,708 KiB | 779.6244507 K; liquid fraction 0 | −1.34106354e-29 J |
| Optimized, same restart and interval | 271.030 s | 781,580 KiB | 779.6244507 K; liquid fraction 0 | −1.34106354e-29 J |

The optimized run is 5.84% faster in this single short measurement. Temperature, enthalpy, liquid fraction and tracker-label arrays at 0.36 ps are bitwise equal to the baseline (maximum absolute difference 0 for each). The optimized 0.36 ps checkpoint also passed an independent restore-only validation. This timing gain is an observed result for one short segment, not a guarantee of the same end-to-end speedup.

The matched solid-state 0.35 ps checkpoints also exercised the 1×/2× comparison pipeline: peak Ag temperature was 809.698 K at 1× and 777.189 K at 2×, while the Ag-only block-averaged temperature RMSE was 1.954 K. Both were fully solid. The peak difference shows why maximum temperature and block-averaged field agreement are reported separately; it is not a melt-window convergence result.

## Implementation change

Profiling showed repeated whole-grid numpy.isclose scans in state_from_enthalpy_density even though the Case A mesh has no cells eligible for the special no-Ag-at-melting-enthalpy branch. The solver now precomputes the candidate indices once and performs the same comparison only at those indices. This changes neither equations nor discretization. The original 0.35 ps checkpoint and the baseline run are retained unchanged.

## 1× comparison reference

The already completed 1× results are read from ../ag70vol_complete_20261028/MELT_FREEZE_METRICS.csv; they were not rerun. That run reported first nonzero liquid at 4.12 ps, first multivoxel cluster at 4.22 ps, and at 4.52 ps a 65-voxel main cluster with physical volume 1864.97 Å³. Its peak was 1264.38 K, Ag liquid volume fraction 1.89234e-4, and relative energy residual −1.61e-12. The reported liquid state became fully solid by 4.548 ps after heat-off at 4.52 ps. These are comparison targets, not evidence that the 2× case will match them.

The 2× onset times, fields, cluster geometry and energy residuals remain pending until the active continuation reaches the corresponding window.

The first long-run segment covered 0.36–0.41 ps in approximately 1,354 s. Its relative energy residual was −3.10e-13, and the saved 0.41 ps checkpoint independently restored with matching temperature, liquid field, tracker time and energy ledger. The live job has continued beyond that restart point.

## Active run and recovery

The long continuation started at the validated optimized_0p35_to0p36/checkpoints/caseA_ag70_2x_final_0p36ps.npz. Each completed heating segment writes a checkpoint, segment CSV and flushed JSONL ledger under heating_0p36_to_pool/. The latest complete validated checkpoint is currently caseA_ag70_2x_heat_0p41ps.npz; before any restart, validate the latest complete checkpoint with continue_ag70_heating.py --restore-only and never restart from a partial temperature log.

Reproduction command from the repository root:

    python -u models/ag_ti3sic2_qsgs_1to1_finer/runs/single_temperature_mesh2_validation_20261027/continue_ag70_heating.py --refinement 2 --restart-from models/ag_ti3sic2_qsgs_1to1_finer/runs/single_temperature_mesh2_validation_20261027/ag70vol_2x_melt_window_20261028/optimized_0p35_to0p36/checkpoints/caseA_ag70_2x_final_0p36ps.npz --max-time-ps 5.0 --max-temperature-K 2300 --pool-followup-ps 0.30 --output-dir models/ag_ti3sic2_qsgs_1to1_finer/runs/single_temperature_mesh2_validation_20261027/ag70vol_2x_melt_window_20261028/heating_0p36_to_pool

## Interpretation limits

The candidate material properties remain incompletely source-verified. Case A omits Ti₃SiC₂ heat storage and conduction and uses insulated exterior boundaries; its cooling is internal energy redistribution, not heat rejection to an environment. Vaporization is absent. A tracked liquid cluster is not by itself proof of a thermally stable melt pool. A clear 1×/2× discrepancy will be reported as unresolved grid sensitivity.
