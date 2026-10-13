# Stage67 step1650 meltpool extraction

**Status: PARTIAL.** The maximum liquid-like Ag cluster has been exported for OVITO. The full 8,728-atom liquid-like Ag selection was not exported because the source package contains only the aggregate count and the 1,032 maximum-cluster IDs; it does not contain per-atom PTM/Q6 classifications. The installed OVITO Basic package has no `ovitos` interpreter, and the available OVITO Python bindings require CPython 3.14, which is not available in this environment.

No LAMMPS run was started. The source dump was read only and was not modified.

## Files

- `meltpool_cluster_step1650.dump` — OVITO-loadable LAMMPS dump of the known 1,032-atom maximum liquid-like Ag cluster.
- `meltpool_cluster_step1650.xyz` — OVITO-loadable extended XYZ of the same cluster, including atom IDs and periodic box metadata.
- `STABLE_MELTPOOL_LIQUID_AG_EXTRACTION_REPORT.md` — counts, provenance, PBC status, and the outstanding full-liquid extraction limitation.
- `input_description.md` — source-frame identity and SHA256; the 31.5 MB original dump is not duplicated in this result folder.
- `SHA256SUMS` — SHA256 hashes of the files shipped here, excluding this manifest itself.

`liquid_Ag_step1650.dump` is **not present**. Creating it would require the per-atom PTM selection for all 8,728 liquid-like Ag atoms; using the 1,032-atom main cluster as a substitute would mislabel a partial subset as the full liquid-like population.
