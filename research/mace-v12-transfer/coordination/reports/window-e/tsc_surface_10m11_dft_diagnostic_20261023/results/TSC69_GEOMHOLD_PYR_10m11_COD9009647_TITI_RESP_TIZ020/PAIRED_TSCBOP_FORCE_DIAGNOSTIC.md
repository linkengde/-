# Paired DFT–TSCBOP force diagnostic: 69-atom Ti/Ti (10-11) slab

## Result

The frozen-geometry DFT single point converged in **57 SCF steps**. A source-ID-matched TSCBOP-only LAMMPS `run 0` was then evaluated on exactly the same geometry. The full-system force-vector RMSE is **2.1356 eV/Å** (component RMSE **1.2330 eV/Å**); the 95th percentile per-atom vector error is **3.3139 eV/Å**, and the maximum is **3.8228 eV/Å** on Ti source atom ID **3**.

| Group | Atoms | Component RMSE (eV/Å) | Vector RMSE (eV/Å) | P95 vector error (eV/Å) | Maximum vector error (eV/Å) | Source ID at maximum |
|---|---:|---:|---:|---:|---:|---:|
| All | 69 | 1.2330 | 2.1356 | 3.3139 | 3.8228 | 3 (Ti) |
| Ti | 33 | 1.4101 | 2.4424 | 3.5561 | 3.8228 | 3 |
| Si | 12 | 0.3057 | 0.5295 | 1.0143 | 1.0143 | 65 |
| C | 24 | 1.2610 | 2.1841 | 3.0749 | 3.2175 | 12 |

The maximum-error Ti atom (ID 3) has nearest C neighbors at 2.0882 Å and nearest Si neighbors at 2.6810 Å; its nearest Ti neighbor is at 2.9696 Å. These are not shorter than the audited minimum distance of the source geometry. The reported DFT maximum force is 2.6229 eV/Å. This is a fixed, as-cleaved surface geometry, so the large forces describe its response at these coordinates and do not indicate that the DFT SCF failed.

## DFT reference and binding

- Structure: `TSC69_GEOMHOLD_PYR_10m11_COD9009647_TITI_RESP_TIZ020`; 69 atoms, C24Si12Ti33, fully periodic; source family `surface_pyramidal_10m11_cod9009647_symm_ti_69atom`, derived from COD 9009647.
- Input SHA-256: `6f7ea4c8584ee36ae4e12488829e8c5ad856a47d4ce1db2772fcb9ff9945cc3d`.
- Geometry SHA-256: `b4e589ceede7089a96c24d4b8bfe260202cfde393c171304095817b56eae599d`.
- Fixed-geometry PBE/PAW, PW500, standard Monkhorst–Pack 6×6×1 without an added Gamma flag or offset; Fermi–Dirac width 0.10 eV; non-spin-polarized; Davidson `niter=2`; MSR1 mixer and SCF energy/density/eigenstate thresholds 1×10⁻⁵; maximum 160 iterations. GPAW 26.7.0, ASE 3.29.0; four MPI ranks and one OpenMP thread.
- Configuration SHA-256: `15ddc8f596281298656a5a21df67b2915f0aebd5ba3bb18316fd2935ee8729c2`.
- Execution-lock SHA-256: `b5a6468db4050a516d2d2ab4d9ba026b9e70076fa682f00eed004739c442f2cc`.
- Runner SHA-256: `b5671637e47a4a179cf17314212fc57a9909a4bc0a65ee5c021fa6d488eaf299`.
- Extrapolated total energy: **−532.1874521237 eV**; force-consistent free energy: **−532.7240450246 eV**. The 69×3 forces and 3×3 stress are finite and present. Input and labelled output match exactly in species/order, Cartesian positions, cell, PBC, and `source_atom_id`.

## TSCBOP-only calculation

- Original potential file: `TSCBOP.source.tersoff`, SHA-256 `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7`.
- LAMMPS 22 Jul 2025 Update 4; `units metal`; `pair_style tersoff`; mapping `Ti Si C`; force units eV/Å; one OpenMP thread; `run 0` only. No minimization and no time integration were run. The potential file was copied unchanged.
- The unusual source-cell frame was rotated only to LAMMPS restricted-triclinic coordinates. Cell metrics and fractional positions were preserved, source atom IDs remained unchanged, and the maximum periodic coordinate round-trip difference was 2.6×10⁻¹⁴ Å. The LAMMPS log confirms 69 atoms and zero steps.
- The absolute LAMMPS potential energy is not compared with the DFT energy; only the same-geometry forces are compared.

## Scope and interpretation

This is a same-parent-family diagnostic, not an independent crystallographic-source validation. The geometry remains `UNASSIGNED_SAME_SOURCE_FACET_DIAGNOSTIC`; it was not assigned to TRAIN, DEV, or TEST and was not used for parameter choice or fitting. The comparison shows sizable TSCBOP force discrepancies in this frozen surface configuration, particularly on Ti and C, but does not by itself establish performance across other surfaces, strains, or thermal states.

This closes the current task. No additional DFT, parameter fitting, MD, or TTM is started automatically. Any later fitting decision still needs independent parent-family DEV/TEST evidence and a separately defined scope.

## Files

- `RESULT.json`, `RUN_STATUS.json`, `labelled.extxyz`, `gpaw_attempt_0000.log`, launcher streams, and their local `SHA256SUMS.txt` are under this task's result directory.
- `tscbop_run0/output/FORCE_COMPARISON.csv` contains one row per source ID, with DFT/TSCBOP force components and vector error.
- `tscbop_run0/output/FORCE_COMPARISON_SUMMARY.json` contains the machine-readable metrics.
- `tscbop_run0/output/RUN0_MANIFEST.json`, `run.in`, `atoms.data`, `coordinate_map.csv`, `forces.dump`, and LAMMPS logs preserve the static comparison provenance.
