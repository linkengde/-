# TSCBOP force-error baseline and root-cause audit

## Auditable reports found

No separately named “Stage E2” root-cause report is present in the verified E commit tree. The closest material is E's `TSCBOP_force_error_environment_analysis.md` and `TSCBOP_literature_applicability_audit.md`, plus A's `candidate02_TSCBOP_bulk_check/{README.md,check.py,report.json,SHA256SUMS.txt}`. The E report's five listed report hashes pass its SHA-256 manifest. The A report and script hashes pass A's own manifest. The TSCBOP file at `research/mace-v12-transfer/coordination/reports/window-b/current_heating_potential_handoff_20261009/potentials/TSCBOP.tersoff` has SHA-256 `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7` and 27 triplet records. It was not modified.

The user supplied DOI `10.1103/PhysRevB.100.214114` and states that the 27 records were numerically compared with Tables V–VII. This audit did not independently retrieve the article. E's literature audit likewise marks external publication verification as unavailable. The parameter file header identifies Plummer and Tucker and Ti₃SiC₂, but file header and numerical table agreement alone do not authenticate the external source.

## Archived numerical baseline

The archived A check compares TSCBOP-only forces with the stored GPAW/PBE force arrays on the same fixed 48-atom cell. Its score is `sqrt(mean_atoms(sum_xyz((F_TSCBOP − F_DFT)^2)))`, in eV/Å. Energy zero points are not matched and no raw absolute energy error is reported.

| DFT mesh on same geometry | All-atom vector force RMSE (eV/Å) | Ti (eV/Å) | C (eV/Å) | Si (eV/Å) | Max per-atom error (eV/Å) |
|---|---:|---:|---:|---:|---:|
| 4×4×2 | 1.996540 | 2.221397 | 2.134622 | 0.003263 | 2.720654 |
| 6×6×2 | 1.987192 | 2.210986 | 2.124643 | 0.003718 | 2.707915 |

The shared input SHA, ordered species/positions/cell, and DFT result manifests were verified in this task. The equality of the three stored baseline geometries (k4, k6, k10) was independently checked from their extxyz contents; only the k4/k6 pair is part of the archived TSCBOP score. The RMS remains near 1.99 eV/Å across those two mesh labels, which are numerical repeats on one geometry, not independent structures.

The archived `check.py` constructs LAMMPS data through ASE with type order Ti/Si/C, runs `units metal`, `pair_style tersoff`, and `pair_coeff * * ... Ti Si C`, sorts returned LAMMPS IDs, checks the 1..48 ID sequence, compares minimum-image positions to ASE coordinates at 1e−10 Å, then scores stored ASE force arrays. This is useful implementation evidence for atom order, mapping, geometry, and eV/Å units. It does not record a separate cell round-trip assertion, does not calculate stress, and does not independently audit every serialization source.

## Reproduction status in this task

This container has neither the ASE package nor the LAMMPS Python module or `lmp` executable. The old `check.py` could not be run here without installing scientific software; per task instruction, no installation was attempted. Therefore:

- archived baseline: verified as a checksummed prior result on two k meshes;
- fresh TSCBOP `run 0` reproduction in this task: not tested;
- unresolved dataset gap: raw 114-atom slab labels are absent from this checkout;
- unresolved physics coverage: no stress/elastic or bulk volume/shear labels are archived for the inspected 48-atom family;
- no confirmed atom-order or force-unit defect was found in the inspected 48-atom archive and checked script, but the independent Stage E2 report itself was not located.

The data/implementation gate is therefore incomplete, and the requested pilot fit is held.
