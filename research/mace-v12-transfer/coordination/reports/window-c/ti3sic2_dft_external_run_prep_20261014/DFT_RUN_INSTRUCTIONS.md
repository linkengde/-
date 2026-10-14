# External DFT Run Instructions — Stage C-D

## Release status

```text
DFT_INPUT_READY=YES
DFT_LABEL_READY=NO
WAITING_FOR_COMPUTE_RESOURCE=YES
SCF_EXECUTION_AUTHORIZED=NO
```

This is a handoff/preparation document only. It does not authorize a calculation. No SCF was run while preparing it. The 15 input geometries remain fixed and unchanged.

## Pinned input source

- Repository: `linkengde/-`
- Input package: `research/mace-v12-transfer/coordination/reports/window-c/ti3sic2_ag_interface_fit_dft_inputs_20261014/`
- Input package source commit: `1fec365ecedb8555ae2c857509d1018953639b86`
- Input manifest SHA-256: `ad15ff9220efd11b6d7288cb45212b853e30b01375a04250f41838954c1b7202`
- Pinned runner: `run_fixed_singlepoint.py`, SHA-256 `d80a83460e35e0d84a964fec763e5afa91808c1543a0978bc83d80d27954e755`
- Current C audit commit: `1c52f2bea37f46bb31aaac20fbe604f4ffbbf297`
- Verify the package `SHA256SUMS.txt` and this handoff's `STRUCTURE_CHECKSUMS.txt` before use.

The 15 geometries comprise six 48-atom Ti₃SiC₂ bulk structures and nine 217-atom Ag/Ti₃SiC₂ Si-termination interface gaps. They are fixed-ion, fixed-cell single-point inputs. Do not regenerate, translate, reorder, relax, or otherwise alter them.

## Software and runtime

The audited environment versions are Python 3.12.14, ASE 3.29.0, GPAW 26.7.0, `gpaw-data` 1.2.1, NumPy 2.5.3, SciPy 1.18.1, and Open MPI 5.0.7. The audited GPAW build reports MPI, ScaLAPACK, and LibXC enabled. Use four MPI ranks and one OpenMP thread per rank (`OMP_NUM_THREADS=1`) for the pinned runner. Confirm the PAW data root through `GPAW_SETUP_PATH`/`gpaw_data.datapath()` and verify the element setup hashes before any authorized run.

Expected PAW setup hashes:

| Setup | SHA-256 |
|---|---|
| `Ag.PBE.gz` | `bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97` |
| `Ti.PBE.gz` | `bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f` |
| `Si.PBE.gz` | `2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0` |
| `C.PBE.gz` | `fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983` |

The package runner is `run_fixed_singlepoint.py`. It validates exact software versions, PAW setup hashes, input/config/runner identities, four MPI ranks, and a separate authorization JSON before calculator construction. All packaged configs have `execution_authorized=false`; no authorization is included in this handoff. Do not bypass or edit the runner gate.

## Frozen DFT method for the 15-point input package

Use the values in each existing JSON config; do not replace them with this prose if a discrepancy is found. Shared settings are:

- Exchange-correlation: PBE; PAW setups from `gpaw-data` 1.2.1.
- GPAW mode: plane-wave, 500 eV cutoff.
- Occupations: Fermi–Dirac, width 0.10 eV.
- SCF criteria: energy, density, and eigenstates each `1e-5` in the config's respective GPAW units; maximum 160 iterations.
- Mixer: Pulay, beta 0.05, `nmaxold=8`, weight 100.
- Fixed ions and fixed cell; no geometry optimization.
- Four MPI ranks, one OpenMP thread per rank.
- Record native total energy and force-consistent free energy separately, complete per-atom force vectors, full stress tensor, SCF convergence/iteration count, software/config/structure hashes, elapsed time, and peak memory.

The six bulk input configs currently specify the default Monkhorst–Pack integer mesh with no explicit Gamma override, at 10×10×4. The nine interface input configs specify a Gamma-centered 2×2×1 mesh. The interface mesh has not been shown converged.

## Two-structure k-point pilot definition

Only after an external resource gate is passed and a separate execution authorization is issued, perform the two requested mesh comparisons on unchanged structures. Within each pair, all settings and the structure SHA must be identical except for the mesh.

| System | Frozen structure | Pilot meshes | Convention and notes |
|---|---|---|---|
| Ti₃SiC₂ bulk | `TSC_VOL_M02`, 48 atoms | 2×2×2 and 4×4×4 | Use the bulk config's default Monkhorst–Pack convention; no explicit Gamma override or extra offset. |
| Ag/Ti₃SiC₂ interface | `AGTSC_SI_GAP_2p4A`, 217 atoms | 1×1×1 and 2×2×1 | Use the interface config's Gamma-centered convention. The 3D slab stress contains vacuum dilution and is diagnostic only. |

These meshes define a proposed pilot comparison, not a modification to the original 15 configs. Preserve original config files. If a separate pilot input is later authorized, it must explicitly bind the same structure/config parent and differ only in k mesh; record its own input hash.

## Comparison and acceptance record

No numerical convergence result is available yet. For each mesh, report:

1. Native total energy per atom and force-consistent free energy per atom as separate quantities; never compare one convention against the other.
2. Maximum atomic force norm in eV/Å, plus the full force array in the result archive.
3. Full stress tensor in eV/Å³. For the interface, label it vacuum-diluted diagnostic data, not a bulk stress target.
4. Difference between the two meshes for each metric, SCF iteration/convergence evidence, memory peak, wall time, and all input/output SHA-256 values.

Do not claim `KPOINT_CONVERGENCE_PASS=YES` until a reviewer has applied a stated numerical error tolerance to the actual paired results. A non-converged SCF, incomplete force array, non-finite output, or missing provenance is not a usable DFT label.

## Archive and stop rules

Keep inputs, outputs, logs, and checkpoints in an isolated run directory. Before launch, verify free memory, cgroup/host limits, scratch capacity, persistent backup capacity, PAW hashes, MPI rank count, and the exact signed authorization. Stop before calculator construction if any check fails. Preserve failed/incomplete logs with an explicit status; do not retry with lower cutoff, changed smearing, changed mixer, relaxed geometry, or altered k-point convention without a separately reviewed change.

The interface pilot is a resource and k-point diagnostic only. The nine-point interface scan and remaining bulk labels are not released by this document. DFT labels remain unavailable until calculations converge and their result archives pass provenance and integrity review.
