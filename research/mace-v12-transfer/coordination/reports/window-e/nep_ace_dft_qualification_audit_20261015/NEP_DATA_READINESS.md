# NEP / GPUMD data readiness

## Format audit

All inspected source `.xyz`/`.extxyz` frames parse with ASE 3.29.0. This proves ASE can read the archived files; it does not prove that the GPUMD/NEP reader accepts their current field names. Existing DFT frames use source-specific keys such as `PW_PBE_energy_eV`, `PW_PBE_forces`, `REF_energy` and related arrays/metadata rather than one verified GPUMD training schema. A future read-only conversion must map scalar energy, species, positions, lattice, PBC and per-atom forces, preserving source IDs and parent-family tags.

Stress/virial is present in only 5 of 71 confirmed energy+force geometries, and its unit/sign convention is not uniform or fully explicit. Do not synthesize missing virials or enable a virial loss until a source-specific conversion rule is checked against the target implementation. Distinguish total energy from force-consistent free energy. Preserve vacuum and PBC exactly; the inventory has `[true,true,true]`, `[true,true,false]` and isolated `[false,false,false]` records.

The requested GPUMD documentation URL was not reachable from this restricted environment (HTTP access denied), and `gpumd` is not installed. Therefore no GPUMD parser or NEP training-file acceptance test was run. This remains unverified.

## Compute environment observed for this audit

The current Codex environment exposes a 4-CPU quota, 32 GiB cgroup memory limit, about 22.4 GB free in `/workspace`, and no visible NVIDIA device or `nvidia-smi`. `gpumd` and `lmp` are absent. These are properties of this audit environment, not claims about every project cloud node. No model inference, training, DFT or MD was run.

## Decision

`NEP_TRAINING_DATA_READY = FAIL`. The archive lacks DFT-labelled Ag liquid and liquid-solid contact families, a clean independent parent-family test, and sufficient consistent stress/virial labels. The parser and GPU requirements have not been operationally verified here. After data additions, first perform schema conversion/round-trip parsing and a small hardware benchmark in an approved GPU environment; do not train yet.
