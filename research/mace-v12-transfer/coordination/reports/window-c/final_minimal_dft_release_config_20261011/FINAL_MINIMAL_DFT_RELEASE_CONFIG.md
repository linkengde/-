# FINAL_MINIMAL_DFT_RELEASE_CONFIG

**Status: `DEFINITION_PARTIALLY_FROZEN_EXECUTION_HOLD`.** This freezes the computable definitions supported by C's reviewed plan. It does not authorize DFT or MACE training.

The E table named in the request, `DFT_PRIORITY_TABLE_20261011.md`, is absent from both the current `cloud-e-b1-dft` tree (head `39734568216e5a3d9853d7d052d28f9799172598`) and the checked `main` tree (head `e91ef26d175072553b409e08f5d9a826929078c7`). Therefore I have not invented its parent IDs, P2 registry, numeric `d0`, or cell dimensions.

## P1 Ti/C

The frozen candidate set is four same-family DEV endpoints: Ti(4f) z-displacements −0.02/+0.02 Å and C(4f) z-displacements −0.02/+0.02 Å. The crystallographic c axis is aligned with z; the parent cell is fixed and only the mapped atom is displaced. Exact parent source and atom IDs remain unresolved: the earlier C plan says the new DEV parent is unselected and excludes simple variants of the already-scored COD 9009647/Kisi 1998 family. Candidate k meshes are 8×8×4 and 10×10×4; the final mesh awaits convergence evidence.

## P2 Ag–Si

The requested family is a new Si-terminated Ti3SiC2(0001)/Ag interface DEV family. `d0` is defined as the gap between the mean z of the outermost Si termination plane and the nearest Ag plane in the verified source geometry. The −0.10/+0.10 Å points are rigid translations of the complete Ag slab relative to the fixed substrate, with lateral registry, internal coordinates, and cell unchanged. The numeric d0, actual registry, parent source, and in-plane cell are unresolved because E's source table and geometry are unavailable. The planned slab requires at least 18 Å z-vacuum and a +4 Å vacuum convergence check. Candidate k meshes are 6×6×1 and 8×8×1; final mesh is pending convergence.

## Shared calculation definition

PBE/PAW plane-wave GPAW, 500 eV cutoff, 0.10 eV Fermi smearing, project SCF thresholds of 1e-5 for energy/density/eigenstates, and a 160-iteration ceiling. Save native energy, force-consistent free energy, complete per-atom forces, and the full stress tensor. Use bulk Ti3SiC2 stress only after convergence; vacuum-containing slab stress is diagnostic and not a training label. P1 and P2 roles are DEV candidates, grouped by parent family. Ag64 is excluded from TRAIN by this task. TEST remains sealed.

See the JSON for structure IDs, unresolved fields, and release gates. Any field marked unresolved must be filled from a verifiable source before an executable input can be frozen.
