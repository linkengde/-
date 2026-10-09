# Ti3SiC2 bulk k-mesh checks

These are fixed-geometry numerical comparisons on the same 48-atom CIF-derived Ti3SiC2 structure with PBE/PW500 settings. Each paired comparison changes only the stated k-point mesh. They do not validate interface forces, independent structures, relaxed cells, or thermal stability.

## Earlier controls

The original 6×6×2 → 8×8×4 diagonal comparison changed both in-plane and c-axis sampling. Its energy differences met the 2 meV/atom budget, while all-atom force-vector RMS was 0.01359 eV/Å and exceeded the 0.01 eV/Å budget. That diagonal comparison could not isolate the mesh direction responsible for the force difference.

## Matched-kz refinements: 8×8 → 10×10

Both new endpoint archives passed their compact checks, and the pair comparisons revalidated SHA manifests, identical input hash, geometry, atom IDs, and method settings other than k points. The per-atom energy differences are signed as 10×10 minus 8×8.

| Fixed kz | Native ΔE (meV/atom) | Free-energy ΔF (meV/atom) | All-atom force-vector RMS (eV/Å) | Max atom-vector difference (eV/Å) | Result |
|---|---:|---:|---:|---:|---|
| 2 | -0.009369 | -0.009753 | 0.006510 | 0.009628 | PASS |
| 4 | -0.035121 | -0.016690 | 0.004618 | 0.007105 | PASS |

Limits: absolute native/free energy difference ≤2 meV/atom; all-atom force-vector RMS ≤0.01 eV/Å.

- A's 10×10×2 run: 37 SCF iterations, four MPI ranks, 3.30 h; archive verification PASS.
- B's 10×10×4 run: 34 SCF iterations, four MPI ranks, 3.26 h; archive verification PASS.

Both fixed-kz refinements meet the current energy and force-vector budgets for this single frozen baseline structure. The earlier 6×6→8×8 force differences were above budget; these 8×8→10×10 results show the increment is smaller at both kz values. They do not establish broad mesh convergence for displaced/strained structures, other geometries, interfaces, or thermally sampled states. Review the full existing mesh matrix before authorizing more calculations.

Machine-readable pair results are `k8x8x2_k10x10x2_comparison.json` and `k8x8x4_k10x10x4_comparison.json`. Recompute either from the repository root with the included `compare_mesh_pair.py`, the respective pair's `--left` and `--right` labels, and an `--output` path.

## Consolidated fixed-geometry matrix

All seven pairwise contrasts were independently recomputed from the published compact archives. Delta means the right-hand mesh minus the left-hand mesh. Each endpoint passed archive verification; the geometry, atom IDs, input hash, and method matched apart from k points.

| Pair | Native ΔE (meV/atom) | Free ΔF (meV/atom) | Force-vector RMS (eV/Å) | Max atom-vector Δ (eV/Å) | Result |
|---|---:|---:|---:|---:|---|
| 4x4x2 → 6x6x2 | 1.735858 | 1.268363 | 0.014596 | 0.020090 | force RMS over budget |
| 6x6x2 → 8x8x2 | -0.163346 | -0.114671 | 0.014436 | 0.022319 | force RMS over budget |
| 8x8x2 → 10x10x2 | -0.009369 | -0.009753 | 0.006510 | 0.009628 | PASS |
| 6x6x2 → 10x10x2 | -0.172715 | -0.124424 | 0.008141 | 0.013019 | PASS |
| 6x6x4 → 8x8x4 | -0.106530 | -0.091330 | 0.010850 | 0.016824 | force RMS over budget |
| 8x8x4 → 10x10x4 | -0.035121 | -0.016690 | 0.004618 | 0.007105 | PASS |
| 6x6x4 → 10x10x4 | -0.141651 | -0.108019 | 0.006258 | 0.009724 | PASS |

The force-vector RMS threshold is 0.01 eV/Å; each energy threshold is an absolute difference of 2 meV/atom. The cumulative 6×6→10×10 comparisons and the final 8×8→10×10 increments pass at both kz=2 and kz=4. Since the 10×10 endpoint is within budget of both 6×6 and 8×8 for this geometry, no further pristine-baseline mesh point is warranted from these results alone. Earlier coarser increments still fail the force budget, so the conclusion is about the 10×10 endpoint on this frozen baseline, not about every lower mesh.

### Next validation step

Use the bulk result only as a numerical reference check. The requested interface-force and separation-force errors remain unmeasured by this matrix. Next, finalize matched slab settings (vacuum size, dipole treatment, in-plane mesh, smearing), then evaluate the existing candidate potential against DFT forces and separation-force projections on interface structures. Keep the gap-scan inputs disabled until those slab settings are fixed; do not infer interface accuracy from the bulk test.

The reproducible matrix is `mesh_matrix_kz2_kz4.json` and `mesh_matrix_kz2_kz4.csv`; regenerate from the repository root with `/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/coordination/reports/window-a/v19_Ti3SiC2_bulk_mesh_comparison/build_mesh_matrix.py`. The script calls the included `compare_mesh_pair.py` for each pair and fails if any archive identity check is not met.
