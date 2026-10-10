# AgSi slab boundary, dipole, and vacuum controls

The original-cell (16 A geometric empty-z) dipole-off and dipole-on controls are complete and archived. Those results are retained as the original-cell comparison.

The next pair was matched to the converged 26 A AgSi slab parent, `AgSi_rigid_Ag_z_parent_k6x6x1_sigma0p10_dipoleXY_v19` (52 atoms; fixed geometry; 6x6x1; PBE/500 eV; sigma 0.10 eV; PBC `[true,true,false]`). The 26 A OFF input is byte-identical to the parent geometry and changes only `poissonsolver` from `{"dipolelayer":"xy"}` to `{}`. It converged in 51 SCF iterations (10,699.7 s); all 18 archive checks passed, and the compact archive was published in commit `7692a6e`. The 36 A ON input expands only the z cell by 10 A and translates all atoms by +5 A; it preserves the slab geometry, IDs, composition, and dipole-ON method. It converged in 54 iterations (13,854.6 s); all 18 archive checks passed, and the compact archive was published in commit `2e37213`.

Run one four-rank GPAW calculation at a time. The 26 A OFF vs 26 A ON pair isolates the dipole-correction effect for that geometry and cell. The 36 A ON vs 26 A ON pair isolates the added vacuum at fixed dipole setting. Compare native energies and free energies separately, per atom, plus same-order force-vector RMS/max. These are fixed-geometry numerical controls, not relaxation, independent validation, training labels, or a claim that either method is converged.

The queue verifies source and input hashes, geometry transform, atom IDs/order, method, owner, MPI/ScaLAPACK, and disk reserve before launch. Compact archives are verified and published; `.gpw` checkpoints remain local and may be removed only after verifier PASS and when needed for disk reserve. `report.json` records the original static-entry audit before runtime integration; this README and the task progress record the subsequent physical run.

## Matched dipole and vacuum comparison

The 26 A dipole-OFF result was compared against the archived 26 A dipole-ON parent. Geometry, cell, atom IDs, method, and ordering match exactly. Turning the correction off changes native energy by −0.000737 meV/atom and free energy by −0.000941 meV/atom; force-vector RMS/max differences are 0.000845/0.002078 eV/Å.

The 36 A dipole-ON result was compared against that same 26 A dipole-ON parent. The registered transformation adds 10 A to the z cell and translates every atom by +5 A; the aligned geometry residual is below 4e-15 A. Native/free energy changes are −0.006762/−0.006551 meV/atom; force-vector RMS/max differences are 0.001056/0.002063 eV/Å.

Both pairs are within the 2 meV/atom energy and 0.01 eV/Å force-vector RMS budgets for this fixed AgSi parent geometry. This supports continuing the AgSi rigid separation scan with 26 A vacuum and dipole correction enabled. It does not establish convergence for other interfaces, positions, or thermal structures. Machine-readable metrics and provenance are in `boundary_comparison.json`.
