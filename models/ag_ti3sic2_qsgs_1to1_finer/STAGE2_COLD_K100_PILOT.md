# Stage 2: K=100 restrained cold transient pilot

Date: 2026-10-24  
Status: completed as a short sensitivity run; **not** a validated equilibrium state and **not** a heating/TTM run.

## Purpose and setup

After the Stage 1 force screen failed, a brief cold transient was run to check whether one strongly restrained configuration could at least be integrated without immediate atom loss or numerical blow-up. This was a diagnostic only. The source structure and potential files were not changed.

- Model: repository QSGS Ag–Ti₃SiC₂ structure, 917,378 atoms; LAMMPS 22 Jul 2025 – Update 4; `units metal`, `boundary p p f`.
- Interactions: the existing Ag EAM, Ti–Si–C TSCBOP/Tersoff, and Ag–Ti/Si/C Morse calibration-v0 terms.
- Fixed layer: atoms initially in z=0–20 Å, held with zero force.
- Thermostat layer: z=20–40 Å, Langevin target 300 K, damping 0.1 ps.
- Mobile region: all atoms outside the fixed layer, integrated with NVE.
- Restraint: `spring/self` with K=100 eV/Å² on the mobile Ti/Si/C skeleton. This is an artificial positional restraint; it does not establish unconstrained skeleton stability.
- Preparation: FIRE minimization capped at 100 iterations with the restraint active; velocities then initialized at 300 K.
- Transient: 250 steps at 0.1 fs (0.025 ps), with no electron system, TTM, or external heat source. The 0.1 fs step was chosen for this short diagnostic and does not validate the previously discussed 0.5 fs production step.

The minimized endpoint before velocity initialization had a maximum Ag force of about 3.69 eV/Å and a net restrained mobile-skeleton force of about 0.208 eV/Å. The former remains above the Stage 1 screening target of 1 eV/Å. Spring forces also mean that restrained skeleton force is not the unrestrained physical force.

## Observed temperatures

The valid run segment reports the following group temperatures. They are transient values, not equilibrium estimates.

| Dynamic time (ps) | Mobile (K) | Ag (K) | Ti/Si/C skeleton (K) | Bath (K) |
|---:|---:|---:|---:|---:|
| 0.005 | 141.5 | 300.9 | 95.1 | 149.3 |
| 0.010 | 152.6 | 306.2 | 107.9 | 154.9 |
| 0.015 | 167.4 | 314.9 | 124.4 | 188.9 |
| 0.020 | 282.7 | 327.9 | 269.5 | 280.4 |
| 0.025 | 156.2 | 346.0 | 101.0 | 180.4 |

The group temperatures vary substantially over only 25 fs. The bath damping time is 0.1 ps, four times the pilot duration, so this run is too short to establish thermal equilibration. The irregular temperature history is a reason to continue diagnosing the cold baseline, not evidence of a material response to heating.

LAMMPS retained all 917,378 atoms through the final step; the run completed without a lost-atom or neighbor-list error. At the endpoint, the reported maximum total per-atom force magnitude was 22.66 eV/Å for the mobile group and 19.13 eV/Å for the mobile skeleton. These values include active spring and Langevin contributions and must not be read as pair-potential residual forces. No claim of force convergence follows from this transient.

## Decision

This pilot only shows that the explicitly restrained model can be integrated for 0.025 ps without immediate loss of atoms. It does **not** pass a cold-equilibrium gate: the duration is too short, group temperatures are not settled, the K=100 positional restraint changes the model response, and the minimized Ag force remains about 3.69 eV/Å. No TTM heat source was started. Do not interpret this as a stable Ag melt-pool result or as validation of the interface potential.

The short run therefore does not authorize a 100 ps heating calculation. The remaining useful work is to keep the no-heat diagnosis explicit and decide whether a constrained sensitivity study is scientifically acceptable despite these limitations; changing the timestep alone cannot remove the force/baseline issue.

## Reproducibility files

The input and large binary outputs are retained in the local run workspace, not committed to the repository:

- Input: `/workspace/runs/ag_ti3sic2_selected/baseline/stage2_cold_k100_sensitivity/in.cold_0p025ps`
- Restart: `/workspace/runs/ag_ti3sic2_selected/baseline/stage2_cold_k100_sensitivity/cold_k100_0p025ps.restart`
- Trajectory: `/workspace/runs/ag_ti3sic2_selected/baseline/stage2_cold_k100_sensitivity/trajectory.dump`
- Raw log: `/workspace/runs/ag_ti3sic2_selected/baseline/stage2_cold_k100_sensitivity/log.cold_0p025ps`

The raw log also contains output from an earlier terminated coordinate-replacement attempt in the same run directory. Its extreme values are stale/invalid and are not used here. The valid run is identified by the `COLD_K100_PROGRESS_0.005PS` through `COLD_K100_0P025PS_DONE_NO_HEAT_NO_TTM` markers and the matching live LAMMPS output. Do not treat the raw log as one clean continuous run.
