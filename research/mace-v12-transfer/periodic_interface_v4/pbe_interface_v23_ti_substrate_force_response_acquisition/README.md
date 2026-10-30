# Ti₃SiC₂ interface Ti force-response acquisition v23

The train49 fit reduced Ag components on some frozen-dev3 interfaces, while the Ti force-vector RMSE increased on all three. This experiment adds six local DFT force labels targeted at the four symmetry-equivalent topmost Ti atoms in the registered AgTi, AgSi and AgC registry-B slab geometries.

## Inputs and geometry audit

For each interface, the four topmost Ti atoms are displaced together by −0.02 Å and +0.02 Å along the surface normal. The other atoms, cell, registry and PBC stay fixed. Inputs contain 52 atoms and no energy or force labels. The parent geometries are the registered v21 interface inputs from the Jeitschko–Nowotny 1967 source family.

The screening parser reads only element identities, positions, cell, PBC and persistent atom IDs. It does not load energies or forces. Against train49, all six inputs have zero exact matches and one intentional near match each: the corresponding parent geometry already in train49. This is expected for a small-displacement force-response acquisition and means these candidates are not independent validation structures. There are no other train49 near hits. Against the frozen dev3 geometry, all six have zero exact or near hits. The official train36/dev3 split is unchanged; test or sealed data was not opened.

The input manifest records the exact parent, dataset, input and source-commit hashes, target atom IDs, displacement, screening results and role. A source geometry repeated in training is used only as the parent for this acquisition; none of the frozen-dev3 structures or labels is copied into the candidate set.

## DFT method and safeguards

Each job is a fixed-geometry PBE single point using GPAW 26.7.0 / ASE 3.29.0, a 500 eV plane-wave cutoff, 6×6×1 k-points, 0.10 eV Fermi–Dirac width, xy dipole correction, four MPI ranks and the registered PAW setup hashes. The native extrapolated energy and finite-width free energy are saved separately; force labels are the GPAW/ASE forces at 0.10 eV width. No relaxation or exact energy-force gradient consistency is claimed.

The environment preflight passed with the 32 GiB cgroup ceiling, 30 GiB stop guard, 0.39 GiB current memory, 28.4 GiB free disk, MPI/ScaLAPACK enabled, verified GPAW/ASE/PAW versions and hashes, and the shared DFT lock available. Jobs run one at a time in AgTi, AgSi, AgC order. At each 20th SCF step, GPAW writes a checkpoint atomically; `checkpoint_keeper.py` retains one hash-verified latest snapshot outside the temporary run directory. If the sequence is interrupted, `recovery_monitor.py` attempts one resume per label, and `resume_sequence.py` validates the checkpoint label and input hash before reuse. If no matching checkpoint exists, the one recovery attempt starts from the registered input. Large `.gpw` files remain local and are not included in result commits. `publish_verified_archives.py` pushes each result archive to `main` only after verification passes.

Run from the repository root:

```bash
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v23_ti_substrate_force_response_acquisition/build_ti_response_inputs.py
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v23_ti_substrate_force_response_acquisition/preflight_environment.py
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v23_ti_substrate_force_response_acquisition/run_sequence.py
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v23_ti_substrate_force_response_acquisition/publish_verified_archives.py
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v23_ti_substrate_force_response_acquisition/checkpoint_keeper.py
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v23_ti_substrate_force_response_acquisition/recovery_monitor.py
```

Run the publisher, checkpoint keeper and recovery monitor in separate processes alongside the sequence. If the automatic resume attempt is exhausted, invoke `resume_sequence.py`; it refuses to run while a DFT process is active and reuses a checkpoint only when its label, input hash and file hash all match.

## Role and next gate

All six results remain `train_acquisition_candidate_pending_audit`. They are not inserted into official training data by this workflow. Once every archive passes `verify_preflight.py`, the follow-up is to build a candidate from the authorized training acquisition pool and rerun the fixed V19 fit against the frozen dev3 only. The acceptance decision remains per-interface force-vector RMSE ≤0.05 eV/Å and separation-force RMSE ≤0.10 eV/Å. A passing fit would still require a separate data-role audit before any formal training integration.

## Status

Inputs, geometry audit and environment preflight are complete, and the sequential DFT run has been launched. `sequence_status.json` is the live machine-readable execution record. Each completed archive is re-verified and pushed to `main` by the publisher. The current run's checkpoint is copied to the local ignored `calculations/latest_checkpoint/` directory once stable.
