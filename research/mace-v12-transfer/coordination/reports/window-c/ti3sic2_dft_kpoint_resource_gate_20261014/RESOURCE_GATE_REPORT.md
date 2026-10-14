# Stage C-D Resource Gate Report

Captured: 2026-10-14 (UTC)  
Repository branch: `cloud-c-e7-tiscbop-readiness-20261014`  
Repository HEAD at audit: `1d428859b313d6a214ea553b5b002fbc5de941d7`  
Input package: `ti3sic2_ag_interface_fit_dft_inputs_20261014`, bound to source package commit `1fec365ecedb8555ae2c857509d1018953639b86`.

## Decision

```text
RESOURCE_GATE_PASS=NO
KPOINT_CONVERGENCE_PASS=NO
DFT_LABEL_GENERATION_READY=NO
SCF_STARTED=0
```

The software stack and MPI support are available, but the current container does not pass the memory safety gate for either requested pilot. No GPAW SCF was launched. The input structures and scientific parameters were not changed.

## Live resource and software evidence

| Check | Observed | Result |
|---|---:|---|
| CPU quota | `cpu.max = 400000 100000` = 4 CPU-equivalents; 5 logical CPUs visible | Available; does not establish exclusive physical cores |
| cgroup memory limit | 34,359,738,368 B (32 GiB) | Hard cap |
| cgroup memory current | 28,327,940,096 B | 6,031,798,272 B (~5.62 GiB) remain below cap |
| cgroup memory peak | 34,359,664,640 B | Historical peak was 73,728 B below cap |
| memory events | `max=15`, `oom=0`, `oom_kill=0`, `oom_group_kill=0` | Prior limit pressure; no recorded OOM kill |
| `/workspace` available | 22,072,561,664 B on overlayfs | Scratch is not a persistent-storage guarantee |
| Free inodes | 2,018,299 | Pass for the small input/report package |
| GPAW | 26.7.0 | Installed |
| ASE | 3.29.0 | Installed |
| Python | 3.12.14 | Installed |
| NumPy / SciPy | 2.5.3 / 1.18.1 | Installed |
| gpaw-data | 1.2.1 | Installed |
| MPI | Open MPI 5.0.7; `gpaw info` reports MPI enabled | Available |
| ScaLAPACK / LibXC | enabled / 5.2.3 | Available |
| GPAW PAW data | Ag, Ti, Si, C setup files present; hashes match the input package records | Pass |
| Active scientific process | No GPAW/DFT process observed during preflight | Pass; no new computation running |

The cgroup high-water mark and 15 `memory.max` events show that this same 32-GiB allocation has approached its hard limit before. Only about 5.62 GiB is currently available under the cgroup cap. The prior E7 resource assessment explicitly disallowed even the 48-atom GPAW points in this environment. The 217-atom interface has no matched successful memory benchmark. Together, these facts do not establish a safe peak-memory reserve for the requested bulk-plus-interface pilot pair.

The host-level `MemAvailable` value is not used to override the container cgroup limit. The workspace overlay free-space figure is also not evidence of persistence or sufficient checkpoint recovery capacity.

## Input integrity and scope

The 15-structure input package manifest was checked with `sha256sum -c`; all listed files passed. The existing package validator reports 15 geometries and no labels or calculations. Relevant pilot candidates are:

| Pilot | Existing structure | Fixed-geometry comparison requested | Status |
|---|---|---|---|
| Ti₃SiC₂ bulk | `TSC_VOL_M02` (48 atoms; existing structure and config hashes retained in its config) | k mesh 2×2×2 versus 4×4×4 | Not run: resource gate failed |
| Ag/Ti₃SiC₂ interface | `AGTSC_SI_GAP_2p4A` (217 atoms; existing structure and config hashes retained in its config) | k mesh 1×1×1 versus 2×2×1 | Not run: resource gate failed |

The interface configuration already states that 2×2×1 convergence is not demonstrated. Its 3D stress is diagnostic for a slab cell containing vacuum and must not be interpreted as a bulk stress target. No input geometry, composition, atom ordering, cell, PBC, or TSCBOP file was edited.

If the resource gate is cleared in a suitably provisioned allocation, each pair must use the exact same frozen structure and all the same DFT settings except the k mesh. Compare like energy conventions: native energy per atom against native energy per atom, and force-consistent free energy per atom against the same quantity. Report maximum atomic force norm and the complete stress tensor with units and the interface vacuum caveat. No comparison values exist from this run.

## Evidence references

- Prior resource basis: `ti3sic2_DFT_extension_20261014/RESOURCE_ESTIMATE.md`.
- Latest structured preflight snapshot: `ti3sic2_ag_interface_dft_labels_20261014/PREFLIGHT_SNAPSHOT.json`.
- Input audit: `ti3sic2_ag_interface_fit_dft_inputs_20261014/DFT_INPUT_REVIEW.md`.
- Input manifest and checksum list: `ti3sic2_ag_interface_fit_dft_inputs_20261014/INPUT_MANIFEST.json` and `SHA256SUMS.txt`.

The previous resource assessment records that a 48-atom baseline at 10×10×4 took 3.26 hours on four MPI ranks; this is not a direct timing estimate for the requested lower meshes. No matched peak-memory result exists for the 217-atom interface, so no safe runtime or memory estimate is asserted for that calculation.
