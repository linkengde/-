# QSGS initial force localization against skeleton morphology

This report correlates the existing full-potential `run 0` force snapshot with each skeleton atom's distance to the QSGS voxel-mask boundary. The snapshot was regenerated from the archived 917,378-atom model with the unchanged repository potentials. No atom positions were advanced; no thermostat, TTM, or heat source was used.

## Finding

High Ti/C forces are present in the voxel-mask interior, not only at thin skeleton surfaces or Ag interfaces. In the `8+ Å` voxel-depth bin, 33,356 Ti atoms have median force 2.67488 eV/Å and 22,193 C atoms have median force 2.04062 eV/Å. These values match the forces on the corresponding Ti and C sites in the pristine parent CIF crystal under the same TSCBOP potential before relaxation. The isolated periodic bulk test archived in [`../ti3sic2_bulk_tscbop/README.md`](../ti3sic2_bulk_tscbop/README.md) showed that coordinate relaxation can remove those bulk forces, although zero-pressure relaxation expands the cell by about 0.5–0.8%.

This is strong evidence that much of the backbone's starting force comes from using the unrelaxed CIF coordinates with the current TSCBOP parameterization, in addition to the separate surface/coordination and Ag–skeleton Morse effects. It does not show that the potential is intrinsically invalid for the crystalline parent phase. Nor does the voxel-depth statistic identify every atomistic surface or bond.

| Element | QSGS skeleton atoms | Median force in voxel depth ≥8 Å | 95th percentile in that bin | Fraction above 2 eV/Å |
|---|---:|---:|---:|---:|
| Ti | 33,356 | 2.67488 eV/Å | 2.67488 eV/Å | 64.7% |
| Si | 11,067 | 0.000049 eV/Å | 0.09399 eV/Å | 0% |
| C | 22,193 | 2.04062 eV/Å | 2.04062 eV/Å | 97.0% |

Across all 704,382 Ti/Si/C atoms, the median force is 1.760 eV/Å, the 95th percentile is 2.709 eV/Å, and 42.1% exceed 2 eV/Å. The backbone-force pattern is therefore not explained by interface forces alone.

## Constrained full-potential relaxation screen

As a follow-up, the 619,397 skeleton atoms above the fixed bottom slab were minimized for 100 FIRE iterations while all 212,996 Ag atoms were held at their original coordinates. The original EAM, TSCBOP, and all three calibration-v0 Morse terms remained active. The minimization hit its iteration cap. Releasing the Ag constraint and evaluating `run 0` gave a maximum force of 3.324 eV/Å on the mobile skeleton and 4.055 eV/Å on Ag. The skeleton moved by at most 0.223 Å (95th percentile 0.088 Å); Ag did not move. The largest Ag force occurs on the same atom ID 717010 that was the initial full-force hotspot.

The minimum Ag–C separation became 2.186 Å, Ag–Ti 2.288 Å, and Ag–Si 2.311 Å. At the repository's 3.25 Å geometric cutoff, the backbone changed from one connected component to three: two C atoms (IDs 175318 and 350658) became isolated, with nearest skeleton distances 3.260 and 3.289 Å. Ag stayed fixed and therefore retained its original connectivity. This constrained relaxation lowers forces for many backbone atoms but does not produce an acceptable continuous, low-force composite state; it is not a cold baseline for heating.

## Method and limits

The saved `qsgs_morphology.npz` mask was tiled three times in x/y before a Euclidean distance transform; z remained open. Atom coordinates were mapped to their mask voxel. 99.804% of Ti/Si/C atoms map to a voxel marked as skeleton; the small mismatch reflects atom/voxel discretization. Depth is a morphology proxy. Per-atom force magnitudes include both TSCBOP and Morse contributions, so this script does not decompose them. The comparison to the parent-crystal forces is diagnostic and does not make a relaxed parent cell directly interchangeable with the composite.

## Reproduction

From `models/ag_ti3sic2_qsgs_1to1_finer/`:

```bash
OMP_NUM_THREADS=4 lmp -sf omp -in potential_validation/qsgs_force_localization/in.force_snapshot -log potential_validation/qsgs_force_localization/log.force_snapshot
gzip -1 potential_validation/qsgs_force_localization/initial_forces.dump
python potential_validation/qsgs_force_localization/analyze_force_surface_depth.py
OMP_NUM_THREADS=4 lmp -sf omp -in potential_validation/qsgs_force_localization/in.fixed_ag_full_potential_skeleton_minimize -log potential_validation/qsgs_force_localization/log.fixed_ag_full_potential_skeleton_minimize
gzip -1 potential_validation/qsgs_force_localization/fixed_ag_skeleton_minimized.dump potential_validation/qsgs_force_localization/fixed_ag_skeleton_minimized.data
python potential_validation/qsgs_force_localization/analyze_fixed_ag_skeleton_minimize.py
```

The generated force snapshot is archived as `initial_forces.dump.gz`; its analysis is in `force_by_skeleton_surface_depth.json`. The constrained-screen endpoint and metrics are archived as `fixed_ag_skeleton_minimized.dump.gz`, `fixed_ag_skeleton_minimize_metrics.json`, and the LAMMPS log. The redundant full LAMMPS data output is omitted; the sorted endpoint dump retains atom IDs, types, coordinates, and forces, and the starting data file remains in the model folder. See `SHA256SUMS.txt` for input and output checksums.
