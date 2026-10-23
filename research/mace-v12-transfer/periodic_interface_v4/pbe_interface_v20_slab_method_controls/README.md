# V20 slab k-mesh and smearing controls

## Purpose and scope

This is a fixed-geometry numerical-method ladder for the existing 52-atom AgSi slab parent. It addresses the unresolved nonperiodic-z slab k-mesh and Fermi-smearing checks before any new model dataset is assembled. It does not resolve the separate need for a source-verified Ti3SiC2 crystallographic parent independent of COD 9009647.

The slab geometry, atom order, IDs, cell, and PBC are byte-matched to the already archived AgSi parent. The parent uses `pbc=[true,true,false]`, a 26 A z cell, and `poissonsolver={"dipolelayer":"xy"}`. Its verified `6x6x1`, sigma `0.10 eV` result is reused as the baseline; it is not recalculated here.

## Registered comparisons

All new jobs use PBE, plane waves at 500 eV, fixed ions, four MPI ranks, identical mixer and SCF thresholds, the same xy dipole correction, and the same input geometry. Only one numerical variable changes in each comparison:

| Label | k-points | Fermi sigma (eV) | Compared with |
|---|---:|---:|---|
| `AgSi_slab_k4x4x1_sigma0p10_v20` | 4x4x1 | 0.10 | Existing 6x6x1, 0.10 baseline |
| `AgSi_slab_k8x8x1_sigma0p10_v20` | 8x8x1 | 0.10 | Existing 6x6x1, 0.10 baseline |
| `AgSi_slab_k6x6x1_sigma0p05_v20` | 6x6x1 | 0.05 | Existing 6x6x1, 0.10 baseline |
| `AgSi_slab_k6x6x1_sigma0p20_v20` | 6x6x1 | 0.20 | Existing 6x6x1, 0.10 baseline |

Only one four-rank GPAW calculation may run at a time. Each exact label is hash-checked, owner-checked, and posted as running before `mpirun`. Large `.gpw` state files remain local and are not archived; compact results are archived after SCF and force verification. A checkpoint may be removed after the archive passes verification.

## Comparison budgets

For the same geometry and ordered atoms, comparisons record native extrapolated energy and free energy separately, force-vector RMS and maximum per-atom vector difference, and the marked Ag-Si pair-force projection. The provisional budgets are 2 meV/atom for each energy convention, 0.01 eV/A for force-vector RMS, and 0.02 eV/A for the marked pair projection. Results characterize only this frozen AgSi parent. They are numerical controls, never V20 training labels or independent validation.

Run one registered label at a time from the repository root:

```bash
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v20_slab_method_controls/run_queue.py AgSi_slab_k4x4x1_sigma0p10_v20
```

Then compare verified outputs with the frozen baseline or another new label:

```bash
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v20_slab_method_controls/compare_results.py AgSi_slab_k4x4x1_sigma0p10_v20 AgSi_rigid_Ag_z_parent_k6x6x1_sigma0p10_dipoleXY_v19
```

The second comparison label is read from the existing, independently archived V19 gap-scan package. New controls do not alter or replace the V19 model or its dataset.
