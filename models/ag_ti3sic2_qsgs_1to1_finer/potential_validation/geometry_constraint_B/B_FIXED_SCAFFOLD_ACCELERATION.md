# Fixed-scaffold Model B acceleration and equivalence check

Date: 2026-10-25. Model B treats Ti₃SiC₂ as an immobile geometric obstacle. The original B input still evaluated TSCBOP forces among the fixed Ti/Si/C atoms on every step, although those forces cannot move the frozen atoms and do not enter the separately defined Ag EAM or WCA cross interaction. This note records a B-only optimization that removes those fixed-fixed calculations.

## Exact model reduction under the rigid-scaffold assumption

`potential.geometry_wca_fixed_skeleton.inc` assigns `zero` to fixed-fixed pairs and retains the same Ag–Ag EAM and all three Ag–skeleton WCA terms from `potential.geometry_wca_candidate.inc`. Since every Ti/Si/C atom is fixed and the B cross terms are explicit pair interactions, omitting the ceramic's internal TSCBOP energy/forces cannot change a mobile Ag force. It does remove the ceramic internal energy constant from the reported potential energy. This equivalence is valid only for the rigid geometric-obstacle model; it does not represent a mobile or thermally responsive Ti₃SiC₂ backbone.

## Results

- At the original B coordinates, the maximum mobile-Ag force is 2.1354904 eV/Å with either include. Comparing force dumps for all 212,800 Ag atoms gave maximum vector difference `1.0×10⁻¹³ eV/Å`; the 95th-percentile difference was zero, and 45 values differed only at roundoff level.
- A 20-step, 0.5 fs no-source NVT comparison from identical minimized Ag coordinates and random seed gave identical endpoint positions, velocities, and forces for all 179,246 mobile Ag atoms. The printed 300 K thermostat temperature, kinetic energy, and maximum Ag force also match exactly. The potential energies differ by a constant because fixed ceramic internal energy is omitted.
- With four OpenMP threads, LAMMPS pair-loop time fell from `39.7894 s` with fixed Ti–Si–C TSCBOP to `1.27883 s` with the rigid-scaffold include, a `31.1×` reduction. Reported per-rank memory fell from `828.1 MB` to `439.1 MB`.
- Repeating the 100-step reference TTM–MD coupling pilot with the fixed-scaffold include produced the same reported temperatures, source energy, and electron-to-lattice transfer as the original B pilot. After correcting the energy-ledger double count, the final residual was `0.101330 eV`, consistent with the no-source NVE drift of `0.101291 eV`. Its end-to-end wall time fell from `279.156 s` to `64.551 s` after this potential reduction and vectorizing the per-atom kinetic-energy handoff.

The model still has the same cold-state drift, arbitrary WCA wall strength, fixed ceramic, Python TTM reference solver, and energy-ledger residual. The acceleration changes throughput; it does not establish a cold baseline or create a melt pool.

## Reproduction

Run from `models/ag_ti3sic2_qsgs_1to1_finer/`:

```bash
OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp \
  -in potential_validation/geometry_constraint_B/in.fixed_skeleton_force_check \
  -log /tmp/log.fixed_skeleton_force_check

OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp \
  -in potential_validation/geometry_constraint_B/in.nvt20_full_skeleton \
  -log /tmp/log.nvt20_full_skeleton

OMP_NUM_THREADS=4 /workspace/.venvs/lammps/bin/lmp -sf omp \
  -in potential_validation/geometry_constraint_B/in.nvt20_fixed_skeleton \
  -log /tmp/log.nvt20_fixed_skeleton
```

The two NVT inputs use the same initial coordinates, random seed, timestep, and 20-step schedule. They write endpoint dumps under `/tmp` for direct comparison. `log.nvt20_full_skeleton` and `log.nvt20_fixed_skeleton` are the captured timing and thermodynamic records.
