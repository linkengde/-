# AgSi slab boundary, dipole, and vacuum controls

The original-cell (16 A geometric empty-z) dipole-off and dipole-on controls are complete and archived. Those results are retained as the original-cell comparison.

The next pair is matched to the converged 26 A AgSi slab parent, `AgSi_rigid_Ag_z_parent_k6x6x1_sigma0p10_dipoleXY_v19` (52 atoms; fixed geometry; 6x6x1; PBE/500 eV; sigma 0.10 eV; PBC `[true,true,false]`). The 26 A OFF input is byte-identical to the parent geometry and changes only `poissonsolver` from `{"dipolelayer":"xy"}` to `{}`. It is enabled and assigned to A. The 36 A ON input expands only the z cell by 10 A and translates all atoms by +5 A; it preserves the slab geometry, IDs, composition, and dipole-ON method. It remains queued until the 26 A OFF archive passes verification.

Run one four-rank GPAW calculation at a time. The 26 A OFF vs 26 A ON pair isolates the dipole-correction effect for that geometry and cell. The 36 A ON vs 26 A ON pair isolates the added vacuum at fixed dipole setting. Compare native energies and free energies separately, per atom, plus same-order force-vector RMS/max. These are fixed-geometry numerical controls, not relaxation, independent validation, training labels, or a claim that either method is converged.

The queue verifies source and input hashes, geometry transform, atom IDs/order, method, owner, MPI/ScaLAPACK, and disk reserve before launch. Compact archives are verified and published; `.gpw` checkpoints remain local and may be removed only after verifier PASS and when needed for disk reserve.
