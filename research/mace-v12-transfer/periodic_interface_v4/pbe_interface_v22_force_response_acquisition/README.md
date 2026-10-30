# Jeitschko 1967 Ag/Ti₃SiC₂ force-response acquisition v22

Train43 lowered the frozen-dev3 aggregate force-vector RMSE from 0.2096 to 0.0600 eV/Å, but all three per-interface force gates remain just above 0.05 eV/Å and AgTi's separation-force error remains 0.1587 eV/Å. This follow-up prepares six fixed-geometry DFT labels that sample the normal force response around the existing AgC, AgSi and AgTi interfaces.

## Registered geometries

Each input rigidly shifts the four Ag atoms together by −0.20 Å or +0.20 Å along the surface normal while holding the source-derived Ti₃SiC₂ slab fixed. The source family and reconstructed 1967 cell are unchanged. Every candidate has 52 atoms, retains the existing PBC and registry, has no energy/force labels, and passed geometry-only screening against train43 and frozen dev3: no exact or near hits were found. The shortest pair distance remains at the 2.118 Å Ti–C backbone distance; the marked Ag-contact distances span 2.527–2.982 Å.

The six inputs are acquisition candidates, not official training assignments. All are from one Jeitschko–Nowotny 1967 source family. The owner-side B review was waived by the user for experiment-only work. The geometry screen reads only atomic identities, coordinates, cell and PBC; it does not inspect labels or test/sealed data.

## Calculation method and safeguards

The fixed-geometry single points use GPAW 26.7.0 / ASE 3.29.0, PBE, a 500 eV plane-wave cutoff, Fermi–Dirac width 0.10 eV, 6×6×1 k points, xy dipole correction, four MPI ranks and the registered PAW setup hashes. Energy is recorded both as native extrapolated energy and finite-width free energy; forces remain the GPAW/ASE forces at 0.10 eV width. No relaxation is performed.

`launch_preflight.py` checks four-rank MPI/ScaLAPACK, the PAW hashes, less than 3 GiB starting memory, at least 5 GiB free disk, the 32 GiB cgroup ceiling with a 30 GiB stop guard, a shared DFT lock and single-job execution. `run_sequence.py` prioritizes AgTi, then AgSi and AgC, runs one job at a time, verifies each archive before proceeding, and never retries automatically. A failed job stops the sequence and preserves its attempt for diagnosis.

Run from the repository root:

```bash
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v22_force_response_acquisition/build_force_response_inputs.py
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v22_force_response_acquisition/run_sequence.py
```

## Current status

All six registered fixed-geometry DFT labels completed on 2026-10-30 at 01:30 UTC (09:30 China Standard Time). Each result passed verify_preflight.py on the completed archive, including the input hash, exact geometry/PBC, method, SCF convergence, finite energies and forces, and internal file checksums.

| Interface | Normal displacement | SCF iterations | GPAW elapsed | Observed peak memory |
| --- | ---: | ---: | ---: | ---: |
| AgTi | −0.20 Å | 54 | 10,514.9 s | 8.79 GiB |
| AgTi | +0.20 Å | 59 | 10,950.3 s | 8.80 GiB |
| AgSi | −0.20 Å | 53 | 11,585.0 s | 8.79 GiB |
| AgSi | +0.20 Å | 54 | 12,598.6 s | 8.74 GiB |
| AgC | −0.20 Å | 61 | 13,054.9 s | 8.72 GiB |
| AgC | +0.20 Å | 62 | 11,399.8 s | 8.75 GiB |

The GPAW summaries record the full calculation elapsed time. For the first AgTi −0.20 Å item, runtime.json additionally records that the supervisor recovered a running job after its launcher exited; its elapsed_s is only the supervisor's recovery interval. The GPAW summary and log contain the full 10,514.9 s runtime. No calculation exceeded the 30 GiB stop guard.

The labels remain train_acquisition_candidate_pending_audit; no result has been inserted into the official training data. The next authorized experiment is a separate train49 fit (train43 plus these six labels), followed by evaluation on frozen dev3 only, using the established train43 recipe. The official train36/dev3 split and sealed test data remain untouched.
