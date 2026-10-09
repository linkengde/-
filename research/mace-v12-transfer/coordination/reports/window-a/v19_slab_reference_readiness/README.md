# AgSi slab reference readiness audit (window A)

**Scope.** Read-only revalidation of existing CIF-derived slab archives and disabled vacuum proposals, to decide the next numerical-reference point. No new DFT was run for this audit. The analysis script rechecks per-archive SHA256 inventories, archive verification status, convergence, source hashes, exact geometry/ID correspondence, and the reported energy/force metrics.

## Existing paired controls

| Comparison | Δ native energy (meV/atom) | Δ free energy (meV/atom) | Force-vector RMS (eV/Å) | Max atom-vector Δ (eV/Å) | Marked-pair projection Δ (eV/Å) | Result |
|---|---:|---:|---:|---:|---:|---|
| 6×6×1 → 8×8×1, σ=0.10, original 16.0003 Å cell, periodic z | −0.181234 | −0.137103 | 0.002767 | 0.004946 | +0.001466 | Pass |
| 16.0003 → 26.0003 Å, 6×6×1, σ=0.10, periodic z, no dipole | −0.014428 | −0.014251 | 0.000923 | 0.002469 | −0.001607 | Pass |
| Dipole off → on at 16.0003 Å, 6×6×1, σ=0.10, z nonperiodic | +0.001740 | +0.001696 | 0.001211 | 0.003321 | +0.001435 | Pass |

Pass budgets used for these *numerical sensitivity controls*: absolute native/free energy ≤2 meV/atom, all-atom force-vector RMS ≤0.01 eV/Å, and absolute marked-pair projection change ≤0.02 eV/Å. The vacuum pair differs by a rigid +5 Å z translation and +10 Å z-cell expansion; the slab, lateral cell, atom order, IDs, markers, PBC and GPAW method otherwise match. The mesh pair changes only k-points. The dipole pair changes only the Poisson dipole correction.

## What these checks support

The published controls provide a provisional starting recipe of 6×6×1 and σ=0.10 and show small separate vacuum-size and dipole sensitivities on this one AgSi geometry. They do not establish mesh or smearing convergence for the nonperiodic-z slab: the mesh and smearing ladders use fully periodic z. The directly relevant parent has **26.0003 Å geometric empty length, PBC=(T,T,F), and dipolelayer=xy**; its current DFT run is a provisional reference point, not a convergence result.

The reference geometry has a maximum force near 0.95 eV/Å. These are fixed-ion numerical comparisons, not relaxed-interface stability tests. They do not measure candidate-potential errors, do not establish independent validation, and do not establish AgC/AgTi accuracy. All three terminations derive from the same CIF family.

## Next registered calculation

Run only `AgSi_rigid_Ag_z_parent_k6x6x1_sigma0p10_dipoleXY_v19`, owned by window A, as the exact 26 Å / dipole-on boundary reference. Keep the other eight rigid-separation records disabled. After archive verification, assess the exact-method force response and decide whether to run a same-vacuum dipole-off contrast or the 36 Å dipole-on control before enabling any ±0.20 Å separation labels. Do not interpret the parent result as a force-error pass for the candidate potential.

The 26 Å dipole-off and 36 Å dipole-on proposals remain geometry-only, unlabelled controls. Before enabling any ±0.20 Å gap point, run a target-slab 4×4×1/6×6×1/8×8×1 mesh ladder and σ=0.05/0.10/0.20 controls at fixed geometry and boundary settings, then compare 26 Å and 36 Å with one fixed dipole setting. The existing 26 Å dipole-OFF and 36 Å dipole-ON proposals change two variables and cannot isolate vacuum sensitivity. A full vacuum-by-dipole matrix is needed to quantify their interaction.

## Reproducibility

Run from the repository root with the project environment:

```bash
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/coordination/reports/window-a/v19_slab_reference_readiness/audit.py
```

`audit.json` contains machine-readable archive checks, exact recorded methods, geometry tests, pair metrics, proposal hashes and limitations. `SHA256SUMS.txt` inventories this report bundle. The audit reads existing archived labels only; it does not read sealed holdouts or run GPAW/MPI.
