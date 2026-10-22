# V19 independent development labels

This package registers one fixed-geometry development structure for each of AgC, AgSi, and AgTi. Roles were frozen before any new DFT label was produced. The structures are derived from the repository's uploaded COD 9009647 Ti3SiC2 CIF, use new lateral Ag registries, and have distinct terminating/contact species. They are not descendants of the three V17/V18 development parents and were not selected from the six V19 residual acquisition structures.

All three structures share the same Ti3SiC2 crystallographic source. The source CIF is internally attributed to COD 9009647 / AMCSD0013252 and Kisi et al. (1998), DOI 10.1016/S0022-3697(98)00226-1; its internal crystallographic tags and exact rebuild are verified in the cited audit. A direct fetch of the COD source previously returned HTTP 403, so the remote bytes/checksum remain unverified. This is a source-recorded proposal based on the uploaded CIF, not a claim of independent online source verification.

The geometry audit reads only element identities, coordinates, cells, PBC, and atom IDs from train36, the six V19 residual inputs, and the three public V17 development parents. It does not inspect energy/force fields or sealed test inputs/labels. All three new structures have distinct composition/cell/PBC signatures from the comparison structures; the published exact/near comparator found no comparable match, and the descriptive species-pair RDF fingerprints have low maximum cosine similarity. AgSi and AgTi candidate pair comparison also fails the near-duplicate threshold by a wide margin. See `geometry_audit.json` for hashes and all checks.

The fixed single-point method is PBE / GPAW PW, 500 eV, 6x6x1, Fermi-Dirac sigma 0.10 eV, pbc [true,true,false], and `poissonsolver={"dipolelayer":"xy"}`. This is one consistent provisional slab reference recipe; it is not a claim of complete slab k-point, smearing, or vacuum convergence. Ions are not relaxed. Native extrapolated energy and free energy are recorded separately.

The three labels are owned by window A and launch-enabled only after input hashes, source lineage, data roles, geometry screening, four-rank MPI/ScaLAPACK availability, and disk reserve are verified. Each converged archive is checked and pushed to `main` before the next label starts. Failed or interrupted calculations remain preserved and are not automatically repeated. Training remains blocked until window B independently audits the three archives and A records `A_APPROVED_V19_DATASET`.

## Registered labels

| Label | Role | Contact | Input |
|---|---|---|---|
| `AgC_v19_independent_dev_01` | independent development | Ag-C | `inputs/AgC_v19_independent_dev_01.extxyz` |
| `AgSi_v19_independent_dev_01` | independent development | Ag-Si | `inputs/AgSi_v19_independent_dev_01.extxyz` |
| `AgTi_v19_independent_dev_01` | independent development | Ag-Ti | `inputs/AgTi_v19_independent_dev_01.extxyz` |

## Reproduce geometry and data-role checks

```bash
/workspace/.venvs/gpaw-mpi/bin/python build_independent_development_inputs.py
```

The builder consumes no sealed paths and does not run DFT. The DFT queue is `run_queue.py LABEL`; it verifies the manifest and owner before launch.
