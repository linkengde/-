# Stage 1: static relaxation diagnostics

Date: 2026-10-24  
Scope: screen the current QSGS geometry and `calibration-v0` potential before any TTM heating. These are independent static minimizations only; no MD steps, thermostat, electron system, or heat source were run. The requested 0.5 fs MD timestep therefore does not enter these results.

## Inputs and reproducibility

The source geometry and potential files were left unchanged. Each test ran from a separate directory under `/workspace/runs/ag_ti3sic2_selected/baseline/stage1_interface_diagnostics/` using LAMMPS 22 Jul 2025 – Update 4, four OpenMP threads, `units metal`, `boundary p p f`.

| Input | SHA-256 |
|---|---|
| `ag_ti3sic2_qsgs.data.gz` | `efc3ce9243d1a1905f04f2fada76d096320761c9a88ff29858365e4f0177675e` |
| `potential.interface_calibration_v0.inc` | `34deced5e08174ba369290b87f5c3bfd2b284b272991adf302b2b0d28c12815f` |
| `Ag_u3.eam` | `cb3d605cd27ec6efb50c4dd7462c2ab3d967e63775dc4400de9134315a44d807` |
| `TSCBOP.tersoff` | `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7` |

The current cross-interface Morse terms remain explicitly marked calibration version 0, not a production potential. The Ag–C Morse equilibrium distance is 2.15256 Å.

## Initial force snapshot

The initial structure contains 917,378 atoms. A `run 0` force dump gives a maximum force-vector magnitude of 3.98780 eV/Å on Ag atom 717010 at `(16.4184, 256.9870, 53.1700) Å`; its force is `(3.1563, −2.1644, −1.1205) eV/Å`. Nearby atoms include C 50510 at 2.3898 Å, C 42643 at 2.5047 Å, C 50504 at 2.5940 Å, and Ag 717000 at 2.8885 Å.

To identify its force source, I compared `run 0` forces with the full potential against a second in-memory `run 0` where only the three Morse `D0` values were set to zero, then matched both dumps by atom ID. For Ag 717010, the full force is `(3.1563, −2.1644, −1.1205) eV/Å`; the EAM contribution is `(-0.5343, -0.5303, 1.0748) eV/Å`; the Morse contribution is `(3.6907, -1.6341, -2.1953) eV/Å`. Their vector magnitudes are 3.9878, 1.3122, and 4.5946 eV/Å, respectively. The cross-interface Morse force is the dominant source at this maximum-force Ag atom, while the Ag–Ag EAM force partly cancels it. The largest Morse-only force belongs to Ag 863396 and is 4.7537 eV/Å; this is a cross-force maximum, not the full-force maximum.

Initial minimum pair distances from the force-dump coordinates are Ag–Ti 2.3323 Å, Ag–Si 2.3320 Å, Ag–C 2.3315 Å, Ti–C 2.0871 Å, Ti–Si 2.6803 Å, and C–C 3.0564 Å. The tracked C pair 497966–498104 initially sits 3.0971 Å apart; both have Ti atom 498108 at about 2.176 Å, while their nearest Ag atoms are about 3.245 Å (C 497966) and 2.519 Å (C 498104). These minima alone do not establish a defect, but they provide reference values for the relaxations below.

## Test A: freeze Ti/Si/C and relax Ag

All Ti, Si, and C atoms were frozen; the bottom 0–20 Å slab was also fixed, leaving 180,363 Ag atoms mobile. CG was capped at 200 iterations (`dmax = 0.02 Å`). The run stopped at the iteration limit.

| Measure | Result |
|---|---:|
| Potential energy | −4,883,783.94 → −4,971,607.20 eV |
| Maximum force-vector magnitude on mobile Ag after minimization | 4.0811 eV/Å |
| Mean / maximum mobile Ag displacement | 0.484 / 2.722 Å |
| Minimum Ag–C distance after minimization | 1.9186 Å |
| Tracked C–C pair 497966–498104 | unchanged at about 3.096 Å because the skeleton was frozen |

Fixing the skeleton did not produce a low-force Ag configuration. The minimum Ag–C separation moved below the Ag–C Morse equilibrium distance, while the maximum Ag force remained above its initial value.

## Test B: harmonic skeleton restraints

Each spring constant was tested from the original geometry in an independent copy. The bottom slab was fixed; all other atoms were minimized for 35 CG iterations (`dmax = 0.02 Å`). The spring was then removed and `run 0` was used to report the physical, unrestrained forces at those coordinates. All three minimizations stopped at their iteration caps.

| Skeleton spring `K` (eV/Å²) | Raw max force on mobile Ag after spring removal (eV/Å) | Raw max force on mobile skeleton after spring removal (eV/Å) | Maximum skeleton displacement (Å) | Minimum Ag–C (Å) | Minimum Ti–C (Å) | Tracked C–C pair 497966–498104 (Å) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 3.833 | 3.111 | 0.534 | 2.064 | 1.859 | 2.860 |
| 10 | 3.799 | 3.252 | 0.278 | 2.065 | 1.926 | 2.960 |
| 100 | 3.704 | 3.492 | 0.035 | 2.047 | 2.050 | 3.076 |

For context, while the spring remained active the reported skeleton force maxima were 2.704, 0.789, and 0.189 eV/Å for K=1, 10, and 100. Those are forces with an external restraint balancing part of the physical force. After removing the restraint, the raw skeleton forces returned to 3.111, 3.252, and 3.492 eV/Å. Thus the stronger spring holds the geometry in place but does not yield a low-force unconstrained state. The K=1 and K=10 cases also shorten the tracked C–C distance.

## Test C: unconstrained FIRE control

Only the bottom slab was fixed. FIRE was capped at 35 iterations (`dmax = 0.01 Å`) and stopped at the iteration limit.

| Measure | Result |
|---|---:|
| Potential energy | −4,883,783.94 → −4,935,424.17 eV |
| Maximum force-vector magnitude on mobile Ag | 3.7949 eV/Å |
| Maximum force-vector magnitude on mobile skeleton | 3.1142 eV/Å |
| Maximum mobile skeleton displacement | 0.357 Å |
| Minimum Ag–C / Ti–C distance | 2.0927 / 1.8692 Å |
| Tracked C–C pair 497966–498104 | 2.9065 Å |

This short FIRE control also did not converge. It follows the same tendency of the tracked C–C pair to shorten. An earlier longer CG audit in `POST_MIN_FORCE_AUDIT.md` reached 2.144 Å for that pair and a maximum force near 14.916 eV/Å; that earlier result is a separate run and is not combined with this test.

## Decision

The current geometry/potential combination does **not** pass the Stage 1 gate for heating. None of the tests reaches the project screening target of less than 1 eV/Å for the maximum force on moving atoms. Ag-only relaxation with a rigid skeleton still leaves a maximum Ag force of 4.081 eV/Å and creates an Ag–C contact at 1.919 Å. Harmonic restraints reduce displacement only while they are applied; raw forces return when removed. The unconstrained FIRE control and prior CG audit both shorten the tracked C–C pair.

Do not start the Ag TTM heat-source pilot or the 100 ps melt-pool run from these minimized coordinates. A timestep of 0.5 fs cannot correct this static force/geometry problem. Keep the original model and potential intact. Under the accepted constraint of retaining the current potential, any next calculation would have to be explicitly labeled as a restraint-dependent sensitivity demonstration, not evidence that an unconstrained composite has a stable Ag melt pool. The Stage 1 acceptance checks should be revisited before heating is authorized.

## Run files

- Initial force snapshot and input: `/workspace/runs/ag_ti3sic2_selected/baseline/stage1_interface_diagnostics/`
- Test A input, log, minimized data, and Ag dump: `testA_ag_only/`
- Test B inputs, logs, and final coordinate dumps: `testB_k1/`, `testB_k10/`, `testB_k100/`
- Test C input, log, and final coordinate dump: `testC_unconstrained_fire/`
- Parsed restraint geometry metrics: `stage1_restraint_geometry_metrics.json`
