# TSCBOP relaxed-parent QSGS candidate: validation status

This is a separate exploratory candidate. It does not replace or modify the accepted QSGS geometry in the parent directory. The purpose was to test whether rebuilding the same bicontinuous 1:1-mass morphology from a Ti₃SiC₂ parent cell relaxed with the repository TSCBOP improves the cold starting state.

## Construction

The parent basis was extracted from the 432-atom Ti₃SiC₂ zero-pressure TSCBOP relaxation documented in `../ti3sic2_bulk_tscbop/README.md`. The relaxed orthorhombic cell is 3.08117938 × 5.33685115 × 17.71017349 Å, about +0.775%, +0.776%, and +0.492% relative to the parent CIF cell. The maximum displacement of a parent basis atom was 0.1061 Å.

`build_candidate.py` uses the same QSGS seed (20261024), 100 × 100 × 50 grid, 26-neighbour growth, 0.0032 core probability, and distance-weighted 0.18 growth probability as the finer geometry. It keeps x/y periodic boundaries, an open z surface with 100 Å vacuum, and a 0–10 Å fixed bottom layer. The final material dimensions are 301.956 × 304.201 × 159.392 Å. The Ag FCC reference lattice parameter is 4.09 Å; the Ti₃SiC₂ c axis is along z. Geometry is deliberately left unrelaxed at construction.

The candidate contains 916,534 atoms: 703,734 Ti₃SiC₂ atoms (117,289 formula units) and 212,800 Ag atoms. Its Ag:Ti₃SiC₂ mass ratio is 0.99999814:1. Independent export/geometry validation finds one connected component per phase, rank-2 x/y periodic winding, and thickness spanning for both phases. This validates the generated geometry only, not its force-field stability.

## Checks performed and outcome

All LAMMPS checks below use the existing Ag EAM + Ti–Si–C TSCBOP + Ag-cross Morse calibration-v0 input. The cross potential remains unvalidated.

| Check | Outcome | Interpretation |
|---|---|---|
| Initial full-potential `run 0` | Maximum force on mobile atoms 3.7801 eV/Å | High residual forces remain in the as-built candidate. |
| Same-coordinate force decomposition | With the three Ag–Ti/Si/C Morse terms diagnostically set to `D0=0`, Ag max force falls from 3.7801 to 1.8506 eV/Å and its p95 falls from 1.7403 to 1.0411 eV/Å. The cross-force component exceeds 1 eV/Å on 27.3% of Ag atoms; 54.1% of Ag atoms lie within 4 Å of the skeleton. | The current cross Morse terms strongly affect the interfacial Ag force tail, but are not the only residual-force source. The no-cross test is a diagnostic decomposition only; `D0=0` is not a proposed production setting. |
| Full-mobile FIRE minimization, first 100 steps | Maximum mobile force 3.6145 eV/Å; one C atom falls outside the 3.25 Å skeleton connectivity cutoff | Iteration cap reached; not converged. |
| Continued full-mobile FIRE, total 200 steps | Maximum mobile force 3.7726 eV/Å; skeleton 8 components (12 atoms outside largest), Ag 5 components (10 outside largest); minimum Ti–C 1.8846 Å and Ag–C 2.0984 Å | Further relaxation worsens local connectivity; reject as a cold baseline. |
| Ag-only FIRE with the entire skeleton fixed, 100 steps | Ag maximum force 3.5640 eV/Å; after releasing the constraint, skeleton maximum force 3.4959 eV/Å | The test held the 0–20 Å bottom slab plus the entire skeleton fixed. It preserves phase connectivity in this short static screen but does not remove large residual forces. |
| Continued Ag-only FIRE, 200 more steps | Ag p95 force drops from 1.640 to 1.513 eV/Å, but maximum rises from 3.564 to 3.845 eV/Å. The endpoint has 58 Ag connected components, with 153 atoms outside the largest component; minimum Ag–C distance is 2.0465 Å. Skeleton maximum force after releasing its constraint is 3.700 eV/Å. | FIRE hit the iteration cap. Lower energy did not produce a uniformly lower force or preserve strict Ag connectivity. |
| Rigid-scaffold, no-source Ag NVT at 300 K, 0.5 fs | 25 fs: Ag temperature 369.75 K; maximum Ag force 3.53 eV/Å | Already shows a strong transient despite the 300 K target. |
| Extended rigid-scaffold, no-source Ag NVT at 300 K, 0.5 fs | Two attempts were manually stopped at 60 fs and 80 fs as the trend worsened. At 80 fs the log reports 1067.28 K and maximum Ag force 11.50 eV/Å. | Not a stable cold baseline. The test had no external heat source or TTM. The NVT temperature rise is an observed failure signal; this short test alone does not identify a unique microscopic cause. |
| No-source Ag NVT after 300 Ag-only FIRE iterations | With the 0–20 Å bottom slab and entire skeleton fixed, temperature rises from 300 to 351.44 K in 25 fs; potential energy drops 1205 eV while kinetic energy rises 1192 eV; maximum Ag force ends at 3.719 eV/Å. The largest Ag component contains 212,532 of 212,800 atoms. | Most of the potential-energy decrease appears as kinetic-energy increase. Ag is still structurally relaxing; this is not externally heated or evidence of melting. |
| Same-coordinate local-force audit at the 300-step endpoint | For mobile Ag, 100,755 atoms are within 4 Å of the skeleton. Their force p95/max is 1.741/3.845 eV/Å; the 8,230 Ag atoms more than 8 Å from the skeleton have p95/max 0.111/0.200 eV/Å. | The high-force Ag tail is concentrated at the Ag–skeleton interface, not in the Ag interior. The zero-cross-Morse calculation is a force decomposition only. |

No TTM run, external heat source, melting run, or cooling/solidification run was performed with this candidate. Do not use its trajectories as evidence of a physical Ag melt pool. The production task remains blocked at the cold-start acceptance step.

The full-mobile FIRE endpoint is not a safe restart: it contains broken phase connectivity and high forces. The 25 fs fixed-scaffold NVT endpoint is diagnostic only; it is not accepted as a production restart. No formal TTM restart was created.

At the 300-step Ag-only endpoint, the maximum mobile Ag force is atom 767369 at (88.689, 96.484, 123.224) Å, about 36.17 Å below the top surface plane. Its nearest skeleton atom is C 198397 at 2.349 Å. At identical coordinates its full-force magnitude is 3.845 eV/Å, versus 0.958 eV/Å with cross Morse removed; the cross-force vector magnitude is 3.962 eV/Å. The largest C force is atom 702729, nearest Ag 706218 at 2.521 Å; its no-cross force is zero at the reported precision, while its full force is 3.662 eV/Å. By comparison, the maximum Ti force is atom 551186 at 3.700 eV/Å, of which 2.966 eV/Å remains in the no-cross diagnostic. This separates a strongly interface-localized Ag/C force contribution from a remaining skeleton-internal force contribution; it does not validate either interaction model.

## Files

- `ag_ti3sic2_qsgs.data.gz`: LAMMPS candidate geometry.
- `geometry_report.json`, `independent_geometry_validation.json`: construction and independent geometry checks.
- `relaxed_parent_basis.json`: parent basis and cell transformation.
- `build_candidate.py`, `derive_relaxed_parent_basis.py`, `validate_geometry.py`: build and validation scripts.
- `in.*`, `log.*`, `*_metrics*.json`: LAMMPS test inputs, logs, and summaries.
- `initial_forces.dump.gz`, `forces_no_cross_morse.dump.gz`, `forces_after300_FULL.dump.gz`, `forces_after300_NO_CROSS.dump.gz`, `mobile_minimized_100.dump.gz`, `mobile_minimized_200.dump.gz`, `ag_only_rigid_skeleton_100.dump.gz`, `ag_only_rigid_skeleton_continue_200.dump.gz`, `rigid_skeleton_Ag_NVT_0p025ps_final.dump.gz`, `rigid_skeleton_Ag_NVT_after300_0p025ps_final.dump.gz`: compressed diagnostic snapshots. `cross_morse_force_decomposition.json`, `high_force_localization_after300.json`, and the rigid-scaffold reports are reproducible with their `analyze_*.py` scripts. No large all-atom trajectory is archived.
- `SHA256SUMS.txt`: checksums for the archived candidate files (excluding the checksum file itself).

The candidate's parent CIF and potential files are shared repository inputs: `research/mace-v12-transfer/periodic_interface_v4/v19_cif_parent_structures/source/Ti3SiC2_COD_9009647.cif`, `potential_validation/ti3sic2_bulk_tscbop/`, and `potential.interface_calibration_v0.inc` in the parent model directory. The potential is unchanged.

## Reproduction

Run from the finer model directory. Required Python packages for building/validation are ASE, NumPy, and SciPy; the local LAMMPS checks used the repository's TSCBOP/EAM/Morse files.

```bash
python potential_validation/relaxed_backbone_candidate/derive_relaxed_parent_basis.py
python potential_validation/relaxed_backbone_candidate/build_candidate.py
python potential_validation/relaxed_backbone_candidate/validate_geometry.py
```

The builder emits raw `.data` and `.extxyz` files. After geometry validation, compress the LAMMPS data file if desired; the archived LAMMPS inputs read `.data.gz`:

```bash
gzip -f potential_validation/relaxed_backbone_candidate/ag_ti3sic2_qsgs.data
OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp -in potential_validation/relaxed_backbone_candidate/in.static_force_check -log potential_validation/relaxed_backbone_candidate/log.static_force_check
OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp -in potential_validation/relaxed_backbone_candidate/in.static_no_cross_morse -log potential_validation/relaxed_backbone_candidate/log.static_no_cross_morse
gzip -f potential_validation/relaxed_backbone_candidate/forces_no_cross_morse.dump
python potential_validation/relaxed_backbone_candidate/analyze_cross_morse_forces.py
OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp -in potential_validation/relaxed_backbone_candidate/in.ag_minimize_rigid_skeleton_continue_200 -log potential_validation/relaxed_backbone_candidate/log.ag_minimize_rigid_skeleton_continue_200
gzip -f potential_validation/relaxed_backbone_candidate/ag_only_rigid_skeleton_continue_200.dump
python potential_validation/relaxed_backbone_candidate/analyze_rigid_screen_continue.py
OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp -in potential_validation/relaxed_backbone_candidate/in.rigid_skeleton_Ag_NVT_after300_0p025ps -log potential_validation/relaxed_backbone_candidate/log.rigid_skeleton_Ag_NVT_after300_0p025ps
gzip -f potential_validation/relaxed_backbone_candidate/rigid_skeleton_Ag_NVT_after300_0p025ps_final.dump
python potential_validation/relaxed_backbone_candidate/analyze_rigid_nvt_after300.py
OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp -in potential_validation/relaxed_backbone_candidate/in.static_forces_after300_full -log potential_validation/relaxed_backbone_candidate/log.static_forces_after300_full
gzip -f potential_validation/relaxed_backbone_candidate/forces_after300_FULL.dump
OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp -in potential_validation/relaxed_backbone_candidate/in.static_forces_after300_no_cross -log potential_validation/relaxed_backbone_candidate/log.static_forces_after300_no_cross
gzip -f potential_validation/relaxed_backbone_candidate/forces_after300_NO_CROSS.dump
python potential_validation/relaxed_backbone_candidate/analyze_high_force_localization.py
```

The second `run 0` input sets the three cross-Morse well depths to zero after including the unchanged calibration-v0 file, solely to subtract force vectors at identical coordinates and estimate the cross-term contribution. Never use that diagnostic input for dynamics.

The archived partial NVT logs end at the stated manual-stop points. The NVT input requests 0.1 ps; reproducing the original stopping point requires stopping that test manually when its trend becomes unacceptable. It is intentionally not presented as a successful equilibration protocol.
