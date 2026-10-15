# Inference cost and environment

## Resource preflight

At execution: Debian 13, Python 3.12.14, 5 CPU cores allocated, cgroup RAM limit 32 GiB with about 28.7 GiB available, no CUDA/GPU, `/workspace` 22 GiB free and `/tmp` 13 GiB free after dependency setup. Inodes had ample headroom. Both models ran CPU-only with one Torch thread. No GPAW, MACE, LAMMPS or system-wide environments were modified. Model weights remained in `/tmp`; no checkpoint was added to Git.

## Isolated inference environments

| Model | Package/runtime | ASE | Torch | NumPy | dtype | Load time | Peak RSS |
|---|---|---|---|---|---|---:|---:|
| SevenNet-MF-0 | SevenNet 0.13.0, official commit `209339d4` | 3.29.0 | 2.5.1+cpu | 2.2.6 | float32 | 5.324 s | 1,146,436 KiB (1.093 GiB) |
| MatterSim-v1.0.0-5M | MatterSim 1.0.0, official v1.0.0 commit `afe58a0` | 3.23.0 | 2.2.0+cpu | 1.26.4 | float32 | 2.661 s | 1,452,564 KiB (1.385 GiB) |

## Per-structure forward-pass times

| Model | 48 bulk | 135 full interface | 21 Ag | 114 TSC | Total of four calls |
|---|---:|---:|---:|---:|---:|
| SevenNet-MF-0 | 0.890 s | 2.574 s | 0.435 s | 3.226 s | 7.125 s |
| MatterSim-v1.0.0-5M | 0.385 s | 1.780 s | 0.045 s | 0.848 s | 3.059 s |

Timings are single local CPU observations, not benchmark guarantees. The inference call times were recorded separately from model loading. Results use different model precision/stack versions; float32-vs-float64 discrepancies are small enough for coarse ranking here but not identical numeric conditions. No GPU result is available.
