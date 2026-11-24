# Ag-W V11 ideal bicontinuous gyroid reconstruction

## Purpose

V11 stops trying to force the damaged QSGS atom coordinates through the old atomic-volume CIC capacity allocator. It builds a new, ideal bicontinuous reference geometry while retaining the main material requirements:

- Ag and W are each continuous through the material region.
- Ag/W mass ratio is tuned to 1:1.
- The material box remains 200 x 200 x 200 A, with a 100 A vacuum region above it (LAMMPS box 200 x 200 x 300 A).
- Ag uses an FCC lattice and W uses a BCC lattice.
- No atoms are placed on top of atoms from the other phase; an explicit narrow interfacial clearance band is included.
- The analytic phase geometry is integrated directly onto the TTM finite-volume mesh.

This is an idealized gyroid reconstruction, not a voxel-identical copy of the original QSGS model. The V7-V10 structures and audit files are untouched.

## Geometry-only prototype results

The local builder prototype produced:

| Quantity | Result |
|---|---:|
| Ag atoms | 291,276 |
| W atoms | 170,906 |
| Total atoms | 462,182 |
| Ag/W mass ratio | 1.000002 |
| Box | 200 x 200 x 300 A |
| Minimum Ag-W distance before EFS relaxation | 2.731 A |
| Ag connected geometry | 100% largest component; spans z |
| W connected geometry | 100% largest component; spans z |
| 10 A grid maximum Ag+W fraction | 1.000 |

The atom count is about 4.2% below the V10 full-structure count. The small reduction is from the explicit interface clearance and mass-ratio tuning, not from arbitrary local deletion. Lattice strains are approximately -0.205% for Ag and +0.462% for W.

## New TTM mapping method

The geometry builder writes Ag, W, and void volume fractions from the same gyroid level set that generated the atoms. Therefore each cell obeys:

f_Ag + f_W + f_void = 1

and f_Ag + f_W <= 1 by construction.

Atom CIC is used only for phase-specific lattice kinetic-energy handoff and Ag source deposition. Each atom's CIC weights are renormalized over cells that contain its own phase. Atomic-volume CIC is no longer used as the finite-volume occupancy estimator. This removes the old local overfill/max-flow bottleneck without moving phase volume tens of angstroms between unrelated cells.

The source remains the exploratory Ag-only Gaussian/exponential profile: I0=5 eV/(ps A2), R=30 A, skin depth=20 A, center=(100,100) A, surface z=200 A. Runtime source power is atom-quadrature based and conserved by phase-aware CIC.

## Files

- scripts/build_agw_gyroid_v11.py: builds the LAMMPS structure, 5/10 A geometry fractions, manifest, report, and SHA-256 list.
- ttm/agw_ttm_geometry_mapping.py: geometry-fraction mapping and phase-aware atom CIC.
- ttm/agw_ttm_md_pilot.py: existing Ag-W TTM pilot, now with optional --geometry-fractions-npz.
- in.v11_static_relax: Dai EFS static minimization input.
- in.v11_cold_baseline: 2 ps 300 K NVT plus 0.05 ps NVE input.
- analyze_v11_outputs.py: static/cold diagnostics and restart-manifest hash update.

## Validation status

The geometry prototype and local numerical checks pass. This is NOT yet an MD or TTM pass. The GitHub Actions workflow builds the data file, runs the EFS minimization, runs the cold NVT/NVE baseline, performs the real LAMMPS mapping-only check, then runs matched source-free and low-energy Gaussian-source TTM software-link trials. Outputs are uploaded as the agw-v11-ideal-gyroid-pipeline workflow artifact.

Even if the TTM energy ledger passes, the result remains an exploratory software/physics baseline. The interface electronic conductance is still provisional, the EFS potential is not a fully validated liquid Ag-W potential, and the frozen geometry fractions are not suitable for long-time melt migration. Dynamic phase remapping is required before interpreting large melt-pool evolution.
