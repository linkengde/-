# B independent review of new-owner mesh results

The B 6x6 pilot converged in 43 SCF iterations, 3723.35 seconds with four MPI ranks. Its compact archive and completion status were published by the registered queue. The local checkpoint is retained and excluded from Git.

`review.py` independently verifies the 4x4/5x5/6x6 archives and checksum inventories, requires identical geometry/cell/PBC/species/IDs/pair marking, and reproduces all five metrics in each published comparison. Run from `/workspace/-`:

```bash
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/coordination/reports/window-b/v19_new_owner_mesh_review/review.py
```

| Meshes | Native energy delta (meV/atom) | Vector difference RMSE (eV/A) | Signed pair delta (eV/A) |
| --- | ---: | ---: | ---: |
| 4/5 | -0.400240983 | 0.021587787 | -0.006197856 |
| 4/6 | 0.657314119 | 0.010490677 | -0.015178681 |
| 5/6 | 1.057555102 | 0.018921803 | -0.008980826 |

All energy and pair gates pass; all vector gates fail against unchanged budgets 2 meV/atom, 0.01 eV/A and 0.02 eV/A. The 4/6 result must not be rounded into a pass. The 5/6 difference exceeds 4/6: monotonic force convergence is not established. These remain numerical-only fixed-ion probes, excluded from training or model-validation truth; SCF convergence does not establish mesh convergence.

A's `new_owner_mesh_diagnosis_and_smearing_plan.json` proposes matched 5x5/6x6 at width 0.05 eV. This is a controlled next experiment, not a demonstrated smearing explanation. Launch remains disabled pending explicit registration and a per-record-width-aware runner/verifier; the existing runner hardcodes 0.1 eV. Preserve historical archives and all other method/input settings. Vacuum/dipole checks remain separate. B starts no unregistered DFT, training, or sealed-label reads.

Verification: three archive verifiers PASS, complete per-archive checksum lists PASS, and all 15 independently computed comparison values reproduce the published results within 1e-12 tolerance. No DFT was rerun for this report.
