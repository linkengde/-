# Geometry-constrained Ag melt-pool model B: finite feasibility trial

This is a separate exploratory branch. It does not modify the original QSGS model or the `Morse calibration-v0` include. It tests the proposed simplification: treat the Ti₃SiC₂ phase as a fixed geometric obstacle and keep Ag–Ag interactions under the repository EAM.

## Model definition

The trial reuses the separate relaxed-parent QSGS geometry in `../relaxed_backbone_candidate/ag_ti3sic2_qsgs.data.gz` (916,534 atoms; Ag:Ti₃SiC₂ mass ratio 0.99999814:1). All Ti/Si/C atoms are fixed during the static and NVT tests, along with all atoms below z=20 Å. Ag remains mobile above z=20 Å.

Ag–Ti/Si/C interactions use a shifted WCA-only Lennard-Jones repulsion. Each zero-force cutoff reuses the corresponding calibration-v0 Morse `r0` (Ag–Ti 2.3371 Å, Ag–Si 2.4041 Å, Ag–C 2.1526 Å); `epsilon=0.5 eV` is a trial obstacle-strength parameter. This removes interfacial attraction and does not represent measured adhesion or physical interfacial heat transfer. It is a geometry-only model, not a replacement production potential.

The wall uses `pair_modify pair lj/cut shift yes`; each cutoff is set to the LJ minimum, so energy and force are zero at the cutoff and repulsive below it. The LAMMPS 2025 executable read and ran this hybrid EAM/Tersoff/LJ setup successfully.

## Finite checks and decision

| Check | Result | Decision |
|---|---|---|
| Initial `run 0` | Maximum mobile-Ag force 2.135 eV/Å; process used about 0.83 GB/rank in this environment. | WCA input loads and evaluates without the Morse cross terms. |
| Ag-only FIRE, 100 iterations; skeleton fixed | Mobile Ag p95/max force 0.615/1.452 eV/Å; p95/max displacement 0.097/0.190 Å. Ag and Ti₃SiC₂ each remain a single connected component. | Better behaved in this bounded static screen; FIRE still hit its iteration cap. |
| WCA contacts at the static endpoint | 166 Ag–Ti and 648 Ag–Si pairs are inside their WCA cutoffs; no Ag–C pair is inside the selected 2.1526 Å cutoff. Maximum single-pair wall force is below 0.70 eV/Å. | Confirms the selected wall is currently weakly active and should not be mistaken for calibrated interface binding. |
| One uninterrupted no-source NVT, 300 K target, 0.5 fs, 100 fs | Ag temperature falls to 183.98 K; potential energy rises 2885 eV, kinetic energy falls 2688 eV, total energy changes 197 eV. p95/max mobile-Ag force is 0.742/2.720 eV/Å; p95/max displacement is 0.304/0.799 Å. Ag has 11 components, with 63 atoms outside the largest component (212,737/212,800 atoms). | No catastrophic no-source heating or force runaway, but the system does not equilibrate at the 300 K target, and strict Ag connectivity degrades slightly. This does **not** pass the cold-baseline criterion. |

The 100 fs test is the bounded B-route acceptance attempt. The earlier two 25 fs segments are retained as diagnostics; the continuous 100 fs result is the decision record. The observed temperature decrease and potential-energy increase indicate continued structural adjustment under the rigid-obstacle model. They are not evidence of melting.

## Conclusion

The geometry-obstacle route reduces the static Ag force tail relative to the Morse-based candidate and avoids the rapid no-source temperature rise observed in the Morse test. It is therefore a useful exploratory model for testing geometric confinement. However, this exact WCA parameterization has not established a stable 300 K baseline, and its interface strength is arbitrary. Do not start TTM heating or claim a physical melt pool from this trial. The model also cannot represent ceramic thermal response or Ag–ceramic adhesion/heat transfer because the skeleton is rigid and cross interactions are purely repulsive.

No local melting, TTM, heat source, automatic shutoff, or cooling/solidification run was performed. The B-route trial is complete at the stated finite endpoint; no additional FIRE/NVT extensions are included.

## Files

- `potential.geometry_wca_candidate.inc`: separate WCA candidate potential. The existing Morse include is untouched.
- `in.static_force_check`, `in.ag_minimize_100`, `in.ag_nvt_0p1ps_continuous`: exact LAMMPS inputs.
- `log.*`, `static_screen_metrics.json`, `nvt_0p1ps_metrics.json`: logs and machine-readable summaries.
- `initial_forces.dump.gz`, `ag_minimized_100.dump.gz`, `ag_nvt_0p1ps_continuous_final.dump.gz`: compressed snapshots.
- `analyze_*.py`: reproducible force, displacement, WCA-contact, and connectivity analyses.
- `SHA256SUMS.txt`: checksums for the archived files other than the checksum list itself.

The QSGS data and Ag/Ti–Si–C intra-phase potentials are shared with the parent validation directory. This B trial uses the separate relaxed-parent geometry candidate, not the canonical 917,378-atom model.

## Reproduction

Run from the finer model directory. Use the repository LAMMPS executable or another compatible build:

```bash
OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp -in potential_validation/geometry_constraint_B/in.static_force_check -log potential_validation/geometry_constraint_B/log.static_force_check
gzip -f potential_validation/geometry_constraint_B/initial_forces.dump
OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp -in potential_validation/geometry_constraint_B/in.ag_minimize_100 -log potential_validation/geometry_constraint_B/log.ag_minimize_100
gzip -f potential_validation/geometry_constraint_B/ag_minimized_100.dump
python potential_validation/geometry_constraint_B/analyze_geometry_constraint_B.py
OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp -in potential_validation/geometry_constraint_B/in.ag_nvt_0p1ps_continuous -log potential_validation/geometry_constraint_B/log.ag_nvt_0p1ps_continuous
gzip -f potential_validation/geometry_constraint_B/ag_nvt_0p1ps_continuous_final.dump
python potential_validation/geometry_constraint_B/analyze_geometry_constraint_B_nvt_0p1ps.py
```
