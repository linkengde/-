# LAMMPS–MACE CPU acceleration and full-model feasibility

## Environment

- Container limit: 4 CPU-equivalents (`cpu.max=400000 100000`), 32 GiB RAM, no CUDA device (`torch 2.5.1+cpu`); OVITO is not installed here.
- `/workspace/.local/bin/lmp` is LAMMPS 4Feb2025 and does not include ML-IAP.
- A separate LAMMPS 30Sep2026 build was produced under `/tmp` with ML-IAP, Python, and double-precision Kokkos Serial/OpenMP. It is isolated from the existing LAMMPS installation and repository model files.
- The frozen MACE V17 weight used for all attempts has SHA256 `559b12c024fe8e2155e26ffc7ca48cf9fa1c2682ea37ff9c5367af14a5cc43d3`; its runtime adapter preserves that weight. Element map checked against the checkpoint: C, Si, Ti, Ag. LAMMPS data atom types are Ti, Si, C, Ag and the `pair_coeff` map was supplied in that order.

## LAMMPS run-0 attempts

The non-Kokkos ML-IAP build loaded the model and cell but failed in MACE's inter-layer communication callback: the standard CPU ML-IAP Python data object has no `forward_exchange` method required by this MACE version. Ignoring exchange would make the periodic graph evaluation wrong, so this path was not patched around or used for dynamics.

The Kokkos-enabled build compiled the dedicated exchange interface and successfully initialized the CPU OpenMP backend, but the first MACE `run 0` exited with `munmap_chunk(): invalid pointer` before emitting energy or forces. A single-thread retry failed the same way. Thus no LAMMPS–MACE force or performance result is available in this environment. The LAMMPS ML-IAP Lennard-Jones smoke test did pass, but that only verifies the LAMMPS package/build, not MACE compatibility.

## Consequence for model size

The full 39,486-atom structure's frozen V17 static evaluation was previously completed in float64, with reported peak RSS about 29.04 GiB. That leaves only about 3 GiB under this container's 32 GiB limit, before velocity arrays, neighbor and communication buffers, trajectory writing, and integration overhead. A full-system MD attempt is therefore not a safe or viable use of this environment. No full-size MD was started.

The 1,926-atom representative slab is the feasible model used for local heating. ASE/MACE calculations are CPU-only; measured resource use fluctuated from about 3 to 18 GiB RSS during the run. Trajectory output is Extended XYZ and can be opened in OVITO on a machine with OVITO installed. No LAMMPS result is represented as completed.

## Status

- `FULL_39486_MD_FEASIBLE=NO` in this 4-CPU/32-GiB/no-GPU container.
- `LAMMPS_MACE_FORCE_CROSSCHECK=FAILED_TO_RUN` due to the integration faults above; no forces were accepted from these attempts.
- `ASE_MACE_SMALL_MODEL_AVAILABLE=YES`.
- These environment limits concern the ability to execute the requested dynamics, not the scientific validity or accuracy of MACE V17.
