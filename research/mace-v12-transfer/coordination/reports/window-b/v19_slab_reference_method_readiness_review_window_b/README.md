# AgSi slab reference-method readiness audit

## Scope

This is a read-only audit of the published AgSi CIF-derived mesh, smearing, and dipole archives; the label-free 26 Å and 36 Å proposals; and A's published Ti3SiC2 bulk mesh matrix. The analysis does not launch DFT/MPI, edit calculation inputs or manifests, enable gap scans, or use sealed labels.

The reproducible analyzer rechecks every selected archive's `SHA256SUMS.txt`, `verification.json`, convergence flag, source-input hash, and exact result/source geometry and atom identity. It independently recomputes signed energy and force differences from the compact archives. All 13 AgSi archives, the proposal bundle, and the Ti3SiC2 report inventory verified; archive and proposal verification records are PASS.

## Matched AgSi archive contrasts

Each delta is right minus left. Energy differences are meV/atom. The force metric is the all-atom force-vector RMS,

`RMS = sqrt(mean_i(||F_right[i] - F_left[i]||_2^2))`.

The componentwise RMS is also recorded separately in `metrics.json`. Reference budgets are absolute native and free-energy deltas ≤2 meV/atom and all-atom force-vector RMS ≤0.01 eV/Å.

### k-point ladder at σ = 0.10 eV

| Pair | Δ native energy | Δ free energy | Force-vector RMS | Assessment |
|---|---:|---:|---:|---|
| Γ → 2×2×1 | 45.191164 | 48.186430 | 1.332792 | Fails energy and force budgets |
| 2×2×1 → 3×3×1 | -6.878073 | -5.858293 | 0.129462 | Fails energy and force budgets |
| 3×3×1 → 4×4×1 | -7.043195 | -5.827788 | 0.074783 | Fails energy and force budgets |
| 4×4×1 → 5×5×1 | -0.400241 | -0.050820 | 0.021588 | Energies pass; force fails |
| 5×5×1 → 6×6×1 | 1.057555 | 0.597101 | 0.018922 | Energies pass; force fails |
| 6×6×1 → 8×8×1 | -0.181234 | -0.137103 | 0.002767 | Passes all three budgets |

All ladder endpoints use the same source hash, atoms, IDs, positions, and cell; only k points change. Their result and source extxyz files report `pbc=[true,true,true]`. Thus the 6×6→8×8 result supports that increment for this fully periodic-z calculation only.

### Smearing at fixed mesh

| Mesh and σ change | Δ native energy | Δ free energy | Force-vector RMS | Assessment |
|---|---:|---:|---:|---|
| 5×5×1: 0.05 → 0.10 eV | 0.195030 | -5.172709 | 0.014224 | Native energy passes; free energy and force fail |
| 5×5×1: 0.10 → 0.20 eV | 0.196070 | -21.549363 | 0.037095 | Native energy passes; free energy and force fail |
| 6×6×1: 0.05 → 0.10 eV | -0.036471 | -5.653965 | 0.015410 | Native energy passes; free energy and force fail |
| 6×6×1: 0.10 → 0.20 eV | -0.588691 | -22.062632 | 0.030258 | Native energy passes; free energy and force fail |

The native-energy changes are small, but the recorded free energies and forces depend materially on the smearing width. These matched width scans also use `pbc=[true,true,true]`; they do not establish width sensitivity for the nonperiodic-z slab.

### Dipole OFF/ON at the same slab geometry

The OFF and ON runs have the same source hash, ordered atoms and IDs, positions, cell, k mesh, smearing, cutoff, XC, and `pbc=[true,true,false]`. The Poisson setting is the only method change. They use a 34.2058765 Å z cell and 16.000267 Å geometric empty length (`cell_z - coordinate_span_z`).

| Pair | Δ native energy | Δ free energy | Force-vector RMS | Maximum atom-vector change |
|---|---:|---:|---:|---:|
| Dipole OFF → ON | 0.001740 | 0.001696 | 0.001211 | 0.003321 |

The contrast is below the stated budgets for this geometry and cell. It directly measures the dipole toggle at the 16 Å geometric empty length; it does not establish the same effect at 26 Å or 36 Å.

## Vacuum proposals and boundary-condition separation

The proposal inventory and verifier both pass. The 26 Å proposal has a 44.2058765 Å z cell and `pbc=[true,true,false]`, with dipole correction OFF. The 36 Å proposal has a 54.2058765 Å z cell and the same PBC, with dipole correction ON. Their coordinate-span complements are 26.000267 Å and 36.000267 Å. These are geometric measures, not electron-density-defined vacuum widths. Neither file contains energy or force labels; both remain unowned and launch-disabled.

The proposed 26 Å OFF → 36 Å ON contrast changes both vacuum and dipole treatment, so it cannot attribute a result to either variable. No energy or force conclusion about vacuum convergence is available yet.

## Boundary conditions control what can be concluded

The mesh and smearing result structures are fully periodic in z even though the cell contains a vacuum-like region and uses one k point along z. The OFF/ON slab archives instead use `pbc=[true,true,false]`. Their atoms, IDs, coordinates, and cell match the mesh structures, but the z boundary condition does not. This makes the mesh/smearing data useful diagnostics for the same geometry, not convergence evidence for the nonperiodic-z slab method.

A's Ti3SiC2 bulk matrix is a separate fixed-geometry bulk result. Its 8×8→10×10 contrasts pass the listed energy and force budgets, but it does not certify an AgSi slab/interface method.

## Recommendation

Use **PBE/PW 500 eV, 6×6×1, σ=0.10 eV, `pbc=[true,true,false]`, and `poissonsolver={"dipolelayer":"xy"}`** as the provisional matched slab reference recipe. Record GPAW 26.7.0 and PAW data 1.2.1 with the registered setup hashes. Use GPAW's native extrapolated energy (`force_consistent=False`) as the primary energy and retain free energy as a separate field.

This is the lowest-cost method already represented by an exact nonperiodic-z slab OFF/ON pair and matches the registered slab proposal family. It is **not certified as mesh- or smearing-converged**, so the gap scan should remain disabled until these matched controls are completed:

1. Compare 4×4×1, 6×6×1, and 8×8×1 on one fixed target slab geometry with the same PBC, vacuum, dipole setting, σ, cutoff, and XC. Require both energy deltas and the all-atom force-vector RMS to meet the stated budgets.
2. At the same target slab geometry and boundary settings, compare σ=0.05, 0.10, and 0.20 eV; check native energy, free energy, and force vectors before fixing a width.
3. Isolate vacuum with matched 26 Å and 36 Å calculations under one dipole setting. To measure dipole effects at both widths, complete the 2×2 vacuum-by-dipole matrix. The current proposals alone do not form a matched pair.
4. After fixing the reference method, assess forces and separation-force projections on the intended gap geometries.

## Reproduce

From the repository root, regenerate `metrics.json` with:

```bash
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/coordination/reports/window-b/v19_slab_reference_method_readiness_review_window_b/audit_slab_reference.py
```

The script uses the existing compact archives and ASE reader; it does not instantiate GPAW or MPI. `SHA256SUMS.txt` covers the README, script, and generated metrics. The metrics preserve each archive's exact path, source/output hashes, method, boundary settings, pair deltas, and individual pass flags.
