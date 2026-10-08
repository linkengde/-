# V19 AgSi matched-mesh comparison: Gamma and 2×2×1

## Result

The Gamma-to-2×2×1 change is far above all three provisional sensitivity budgets. This means Gamma is not interchangeable with 2×2×1 for this fixed geometry under the recorded settings. It does **not** establish that 2×2×1 is converged; A's registered 3×3×1 calculation is still needed for that comparison.

| Observable (2×2×1 minus Gamma) | Change | Provisional budget | Result |
|---|---:|---:|---|
| Native energy per atom | +45.1912 meV/atom | 2 meV/atom | Exceeds |
| Free energy per atom | +48.1864 meV/atom | No separate gate; recorded for context | — |
| Force-vector difference RMSE | 1.33279 eV/Å | 0.01 eV/Å | Exceeds |
| Maximum per-atom force-difference norm | 2.64058 eV/Å | Descriptive | — |
| Marked Ag–Si signed projection change | +1.16883 eV/Å | 0.02 eV/Å | Exceeds |

## Verification and interpretation

Both compact archives pass their `verification.json` checks and their `SHA256SUMS.txt` inventories. Each uses the same exact CIF-derived input (`8fe44bee35e17ebfc433e059925a7a85f0dbb4e775090a75ae1883b3eaa6c36c`), 52 atoms, fixed positions/cell/PBC, PBE PW-500 eV, Fermi smearing 0.1 eV, the same SCF/mixer recipe and four MPI ranks. Gamma converged in 51 iterations; 2×2×1 converged in 41. The archived native energies are −392.1041689830 and −389.7542284612 eV/cell; free energies are −392.7609172268 and −390.2552228713 eV/cell.

The vector metric is the per-atom RMSE of the full force-vector difference. The marked-pair value is the 2×2-minus-Gamma change in `(F_Si − F_Ag)` projected along the minimum-image Ag-to-Si vector. The signs therefore indicate the signed change with increasing mesh density, not an absolute error against a converged target.

The differences are large enough that future labels must use a frozen common k-point convention. Do not choose a mesh as truth from this pair alone, and do not mix Gamma and denser-mesh results into training or development labels. The comparison is a fixed-ion numerical probe on an asymmetric, unrelaxed slab; it gives no equilibrium, model-accuracy, or dataset-entry result. Vacuum/dipole sensitivity remains separate. A's 3×3×1 run is active; B did not start another calculation.

The reproducible comparison result is copied here as `comparison.json`; the canonical pilot-package output is `periodic_interface_v4/pbe_interface_v19_cif_convergence_pilot/mesh_comparison_AgSi_COD9009647_pilot_gamma__AgSi_COD9009647_pilot_k2x2.json`.
