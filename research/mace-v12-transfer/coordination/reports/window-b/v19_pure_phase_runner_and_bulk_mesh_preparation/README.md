# Disabled generic pure-phase DFT and bulk mesh draft

Reuses all14 previously verified pure-phase control inputs byte-for-byte. No source inputs regenerated, geometry audit repeated, DFT/mpirun or model operation executed. Twenty proposed records cover four baseline meshes per phase and the six displacement/volume/shear controls per phase at one proposed mesh. All launches false, owners null, blocked_pending_A_registration.

## Generic entry

entry.py binds exact input/source hashes, phase stoichiometry, ordered symbols/local_new_id, TTT PBC, empty Poisson, fixed proposed PBE/PW500/Fermi0.10, registered method/mesh and energy convention. Ag-only and Ti3SiC2-only are supported; no central_pair or all-four-species assumption. Four PBE setup hashes pin available reference setup files; only setups for present elements are used. Metadata pins GPAW26.7.0/ASE3.29.0/gpaw-data1.2.1, SCF/mixer/maxiter, parallelism and energy keys.

Runner preserves reviewed MPI-before-ASE collective-read ordering and full-allocation CPU checks under one-core rank binding. It requires four MPI ranks, ScaLAPACK, exact input path, no existing run,3GiB startup and1.5GiB SCF reserve, one thread/rank. Queue uses an exclusive DFT lock and rejects active GPAW/existing outputs; publishes stage intent before launching one registered job. Native extrapolated energy and free energy are separate; finite-width forces remain consistent with free energy, not automatically the extrapolated energy derivative. No physical stability/zero-temperature conclusion follows from SCF success.

Output uses17-digit positions/cell/energy with actual PBC and local IDs, no marked pair. Verifier requires complete immutable compact inventory, exact geometry/PBC/order/IDs/formula/source/method/Poisson, finite Nx3 forces and energies, MPI4, consistent convergence/iteration/log evidence and force summaries. state.gpw stays local. Queue verifies twice then uses the repository publication helper; A must register label_result_directories to the actual new production entry before enabling. Keep draft input paths/root inference correct when relocating: current ROOT parents5 resolves this report checkout; production root should be explicitly reviewed. Never enable only a manifest without registered owner/jobs and refreshed inventory.

## Phase-specific three-dimensional ladder

For actual cell matrix A, rows B=2pi*inverse(A).T are reciprocal vectors (Angstrom inverse). Candidate n_i=ceil(|B_i|/target), rounded upward to even integers, keeps |B_i|/n_i<=target and Monkhorst-Pack parity fixed. Targets0.40,0.25,0.16,0.12 inverse A; duplicate meshes removed.

| Phase/input | Mesh gradient | Proposed perturbed-control mesh |
| --- | --- | --- |
| Ag conventional2x2x2,32atoms | 2x2x2 → 4x4x4 → 6x6x6 → 8x8x8 | 4x4x4 |
| Ti3SiC2 CIF2x2x1,48atoms | 4x4x2 → 6x6x2 → 8x8x4 → 10x10x4 | 6x6x2 |

The elongated Ti3SiC2 cell needs a smaller z count at comparable reciprocal spacing; neither uses slab6x6x1. mesh_plan.json records exact reciprocal vectors, spacings, raw kpoint upper bounds and raw-work ratios. These ratios are not measured CPU time: symmetry, band/plane-wave counts, SCF and checkpoint resources remain unknown. Before costly launches A should inspect exact irreducible kpoints and memory/disk, use one four-rank job per machine and preserve outputs.

Compare native/free energy deltas per atom, all-atom vector RMS/max differences, unchanged proposed2meV/atom and0.01eV/A vector budgets; no interface-pair gate applies. Symmetric equilibrium bulk forces can be zero at every mesh, so their pass cannot validate forces on displaced/strained controls. Require mesh comparisons on those controls plus separate width/cutoff checks. Stress is optional, not currently exported: if implemented declare ASE eV/A^3, tension-positive, Voigt xx yy zz yz xz xy and verify units/sign/conventions before comparison. Do not pretend unrun stress/energy/force checks passed.

The ASE Ag4.09 starting lattice is not optimized/externally verified experiment; COD input has uploaded/internal provenance. Pure-phase perturbations share parents and are numerical controls, not independent interface validation or automatic training labels. Common physical target remains A's decision. No collision input, thermostat, MD or production changes.

## Static validation

Temporary fixtures are explicitly NONPHYSICAL and are not archived as DFT. Tests cover disabled runner/queue, enabled-but-unassigned owner, wrong composition/IDs/input hash/PBC/Poisson, changed result cell/species/method, nonfinite forces, wrong convention and false convergence; recomputed fixture checksums cannot hide metadata changes. Valid fixture acceptance tests schema only, not physical readiness.

From this directory, use `/workspace/.venvs/gpaw-mpi/bin/python test_rejections.py`. Regenerate entry SHA256 after generated validation/report changes. Runtime calculations/ and caches are excluded from entry inventory but result inventory covers all compact result files. No runtime DFT was tested or authorized.
