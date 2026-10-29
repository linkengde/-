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

The geometry and environment preflights passed on 2026-10-29. Two labels have completed and passed archive verification: AgTi −0.20 Å (54 SCF iterations, 10,514.9 s, peak 8.79 GiB) and AgTi +0.20 Å (59 SCF iterations, 10,950.3 s, peak 8.80 GiB). The current item is `AgSi_JN1967_registryB_gap_m0p20_sigma0p10_k6x6x1_v22`, at SCF iteration 20 at the latest status check. The original interactive launcher exited after its execution session detached while the MPI calculation continued; the recovery supervisor now watches the active `.running_<label>` attempt, enforces the same 30 GiB stop guard, archives and verifies it, then continues serially. Supervisor and per-job output are redirected to durable files to avoid tying calculation lifetime to an interactive terminal. Three labels remain after the current item. Prior same-protocol 52-atom interface calculations in v21 took 12,621–22,563 seconds apiece, so the remaining work will take substantial wall time. Automatic retries remain disabled.

After all six results pass their per-result archive verifiers, the next dataset will be a separate train49 experiment (train43 plus six labels), evaluated only on frozen dev3. Training integration into official train36/dev3 remains disabled. This panel does not establish geometry-specific k-point convergence for all six perturbed slabs or certify production/MD use.
