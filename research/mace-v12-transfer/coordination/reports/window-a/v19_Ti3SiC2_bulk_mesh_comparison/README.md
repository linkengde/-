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
