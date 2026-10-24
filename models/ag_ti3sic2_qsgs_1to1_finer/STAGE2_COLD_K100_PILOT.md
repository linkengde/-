# Stage 2: K=100 restrained cold transient pilot

Date: 2026-10-24  
Status: completed as a short sensitivity run; **not** a validated equilibrium state and **not** a heating/TTM run. This revision supersedes the earlier file references below, which pointed at stale, invalid outputs.

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

## Independent cold continuation

To check whether the short transient could be extended from its actual end state, the valid final dump was mapped by atom ID onto a fresh read of the original data. This restored the original fixed/bath/Newton membership (117,618 / 119,063 / 680,697 atoms). The Ag and moving-skeleton temperatures at step 0 matched the prior endpoint (345.96 and 100.95 K), confirming that velocities and groups were imported consistently. The K=100 eV/Å² spring was then re-anchored at those positions, the Langevin bath was restarted with a new random seed, and another 250 steps (0.025 ps at 0.1 fs) were run. This is an independent restrained continuation, not a mathematically continuous restart because the spring reference and random-number stream were reset.

| Continuation time (ps) | Mobile (K) | Ag (K) | Ti/Si/C skeleton (K) | Bath (K) |
|---:|---:|---:|---:|---:|
| 0.000 | 156.2 | 346.0 | 101.0 | 180.4 |
| 0.005 | 150.7 | 369.8 | 86.9 | 159.3 |
| 0.010 | 133.6 | 400.4 | 55.9 | 151.0 |
| 0.015 | 177.0 | 438.3 | 100.9 | 190.2 |
| 0.020 | 190.6 | 483.8 | 105.3 | 216.8 |
| 0.025 | 195.7 | 536.9 | 96.3 | 206.8 |

The continuation retained all 917,378 atoms and completed with zero dangerous neighbor builds. Ag warmed steadily to 536.9 K, while the skeleton reached 55.9 K at 0.010 ps and ended at 96.3 K. This is further evidence that the model has not reached a common 300 K cold baseline; with the current non-equilibrium state and artificial restraint, this trend cannot be interpreted as a physical melt process.

## Checkpoint and trajectory integrity correction

The first run directory, `/workspace/runs/ag_ti3sic2_selected/baseline/stage2_cold_k100_sensitivity/`, contains mixed logs and stale output from the terminated coordinate-replacement attempt. In particular, its file named `cold_k100_0p025ps.restart` reads as step 450 with group temperatures around 4.0×10⁷ K, while the purported trajectory has only steps 100, 150, and 200. **Neither that restart nor that trajectory is a valid record of the successful 0.025 ps run; do not use them.** The previous revision incorrectly listed them as reproducibility outputs.

The clean rerun and continuation are in a separate directory. For the clean pilot, `final_valid.dump` and `cold_k100_valid_0p025ps.restart` both contain step 350 and 917,378 atoms. Re-exporting the restart and comparing by atom ID gave exact velocities and z coordinates; periodic minimum-image x/y coordinate differences were at most 0.0005 Å, due to text formatting/wrapping. For the continuation, the final dump and re-exported restart both contain step 250 and 917,378 atoms, with the same round-trip agreement. In `extension_0p025ps/`, use only the `*_corrected` input/log/dump/restart outputs listed below; earlier attempts there were interrupted or used the wrong temperature-reporting groups and are not valid continuation results.

## Decision

This pilot only shows that the explicitly restrained model can be integrated for 0.025 ps without immediate loss of atoms. It does **not** pass a cold-equilibrium gate: the duration is too short, group temperatures are not settled, the K=100 positional restraint changes the model response, and the minimized Ag force remains about 3.69 eV/Å. No TTM heat source was started. Do not interpret this as a stable Ag melt-pool result or as validation of the interface potential.

The short pilot and its continuation do not authorize a 100 ps heating calculation. The force gate still fails, and the cold temperatures remain far from equilibrium even under K=100 restraint. No TTM heat source was started. Changing the timestep alone cannot remove the force/baseline issue.

## Reproducibility files

The input and large binary outputs are retained in the local run workspace, not committed to the repository:

- Clean pilot: `/workspace/runs/ag_ti3sic2_selected/baseline/stage2_cold_k100_repro_20261024/` (`in.cold_0p025ps_repro`, `log.cold_0p025ps_repro`, `trajectory_valid.dump`, `final_valid.dump`, `cold_k100_valid_0p025ps.restart`).
- Independent continuation: `/workspace/runs/ag_ti3sic2_selected/baseline/stage2_cold_k100_repro_20261024/extension_0p025ps/` (`in.extension_0p025ps_corrected`, `log.extension_corrected`, `trajectory_extension_corrected.dump`, `final_extension_corrected.dump`, `cold_k100_extension_corrected.restart`).
- Restart round-trip dumps: `restart_roundtrip.dump` and `extension_0p025ps/restart_extension_roundtrip.dump` in the clean pilot directory.
