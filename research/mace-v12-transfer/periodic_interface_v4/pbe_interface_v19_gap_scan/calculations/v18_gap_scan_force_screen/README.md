# V18 force screen on the AgSi rigid-layer gap scan

## Result

The fixed V18 clean-core epoch-79 model fails both provisional force gates on all three archive-verified configurations. This is a same-CIF-family diagnostic, not independent validation or a production qualification.

| Configuration | All-atom force-vector RMSE (eV/Å), gate ≤0.05 | Marked Ag–Si separation-force absolute error (eV/Å), gate ≤0.10 | Result |
|---|---:|---:|---|
| Ag layer −0.20 Å | 2.63466 | 1.89572 | Fail both |
| Parent | 2.63378 | 1.71651 | Fail both |
| Ag layer +0.20 Å | 2.63371 | 1.27048 | Fail both |

The marked-pair metric preserves the project's existing definition, `(F_Si − F_Ag) · unit(Ag→Si)`. The report also records half this value for an equal-and-opposite relative-coordinate convention, and the total +z force on the entire Ag layer. These are different generalized coordinates and must not be conflated.

For the parent, per-species force-vector RMSE is 0.585 eV/Å for Ag, 0.839 for Si, 3.123 for Ti, and 2.735 for C. Distance localization gives 0.929 eV/Å for atoms within 3 Å of any Ag, versus 2.751 eV/Å for atoms more than 7 Å away. The large global error is therefore mainly in the Ti/C framework away from the Ag layer, consistent across the three shifts. Near-interface errors are lower but remain substantial; the Ag–Si pair separation force also misses its gate by 1.27–1.90 eV/Å.

The central finite-difference generalized force for rigid +z translation of the Ag layer is −2.471 eV/Å from DFT free energies and −4.643 eV/Å from MACE, a model-minus-DFT difference of −2.172 eV/Å. The DFT native-energy derivative is reported separately at −2.493 eV/Å. The free-energy derivative is used for the force-consistent comparison.

## Calculation and scope

- All three geometries are fixed-ion structures from one uploaded CIF family; the ±0.20 Å changes rigidly translate the Ag layer. GPAW SCF converged for all points. They are not ionic relaxations.
- The DFT calculations use PBE, 500 eV plane-wave cutoff, 6×6×1 k-points, 0.10 eV smearing, `pbc=[true,true,false]`, `dipolelayer='xy'`, and four MPI ranks.
- Each DFT archive passed all 18 checks in `verify_result.py`. The +0.20 Å point converged in 68 SCF iterations in about 4 h 23 min. Its 2.8 GB local checkpoint was removed only after the archived result passed verification; no checkpoint is in the repository.
- The candidate is the pinned V18 clean-core epoch-79 model. Lineage hashes were rechecked. Inference used Python 3.12.14, PyTorch 2.5.1+cpu, mace-torch 0.3.16, and ASE 3.29.0 on CPU.
- These three closely related geometries are not independent holdouts. The Ag–Si local contact is in the neighborhood of V18 training coverage. Development/test structures and sealed labels were not read, and no model training or tuning was performed.

## Files and reproduction

- `metrics.json`: archive hashes, model lineage, per-configuration energies and force metrics, gates, species and distance localization, and finite-difference force comparison.
- `per_atom_force_errors.csv`: DFT/MACE force vectors and errors for every atom in all three configurations, including position and distance from the Ag layer.
- `analyze_gap_scan.py`: reproduction script; it re-verifies all three DFT archives and the pinned V18 lineage before inference.
- `SHA256SUMS.txt`: checksums for this report package.

From the repository root, reproduce with:

```bash
/workspace/.venvs/mace-v12/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_gap_scan/calculations/v18_gap_scan_force_screen/analyze_gap_scan.py
```
