# Ti₃SiC₂ interface carbon force-response acquisition v24

The fixed V19 seed-45 train55 fit improved aggregate interface force error over train49, but all three per-interface force-vector RMSE values remain above the 0.05 eV/Å gate. This acquisition adds six local force labels targeting the carbon response nearest the marked Ag contact pair, complementing v23's topmost-Ti response labels.

## Inputs and geometry audit

For each registered AgTi, AgSi and AgC registry-B slab, the four C atoms in the contact-nearest symmetry orbit are displaced together by −0.02 Å and +0.02 Å along the surface normal. The selected orbit is the four-atom C z-plane with the smallest mean distance to the marked central contact pair. All other atom positions, the cell, registry and PBC stay fixed. Each input has 52 atoms and contains no DFT energy or force labels.

Parent coordinates come from the registered v22 interface acquisition inputs at commit `5d1e734c16d17dc9b460e3184eab4670ccc2739c`. Candidate screening uses the train55 candidate dataset and frozen dev3 geometry at commit `6b8fe5520679808cd98506eb08b8ea72cd71f61a`. The geometry-only parser reads element identities, positions, cell, PBC and persistent atom IDs; it does not read energies or forces. Every candidate has zero exact training matches and three intentional near matches: its parent structure and the same-interface v23 ±0.02 Å topmost-Ti response pair. There are no other train55 near matches and no exact or near matches against frozen dev3. The paired C displacements are themselves intentional response probes. This means these labels are acquisition data, not independent validation structures. The official split is unchanged, and no test or sealed data was opened.

## DFT method and safeguards

Each job is a fixed-geometry PBE single point using GPAW 26.7.0 / ASE 3.29.0, a 500 eV plane-wave cutoff, 6×6×1 k-points, 0.10 eV Fermi–Dirac width, xy dipole correction, four MPI ranks and the registered PAW setup hashes. The native extrapolated energy and finite-width free energy are saved separately; force labels are GPAW/ASE forces at 0.10 eV width. No relaxation or exact energy-force gradient consistency is claimed.

The geometry and environment preflight passed: 32 GiB cgroup ceiling, 30 GiB stop guard, memory below the 3 GiB start limit, at least 5 GiB free disk, MPI/ScaLAPACK enabled, verified GPAW/ASE/PAW versions and hashes, and the shared DFT lock available. Jobs run one at a time in AgTi, AgSi, AgC order. At each 20th SCF step, GPAW writes a checkpoint atomically. `checkpoint_keeper.py` retains one hash-verified latest snapshot outside the temporary run directory. `recovery_monitor.py` allows one automatic resume per interrupted label; `resume_sequence.py` validates checkpoint label and input hash before reuse. Large `.gpw` files remain local and are excluded from result commits. `publish_verified_archives.py` pushes each result archive to `main` only after `verify_preflight.py` passes.

Run from the repository root:

```bash
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v24_carbon_force_response_acquisition/preflight_environment.py
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v24_carbon_force_response_acquisition/run_sequence.py
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v24_carbon_force_response_acquisition/publish_verified_archives.py
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v24_carbon_force_response_acquisition/checkpoint_keeper.py
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v24_carbon_force_response_acquisition/recovery_monitor.py
```

Run the publisher, checkpoint keeper and recovery monitor in separate processes alongside the sequence. If the automatic resume attempt is exhausted, invoke `resume_sequence.py`; it refuses to run while a DFT process is active and only reuses a matching verified checkpoint.

## Role and next gate

All six configurations remain `train_acquisition_candidate_pending_audit`; the workflow does not add them to official training data. After all archives pass verification, the next step is to construct a candidate from the authorized acquisition pool and rerun the fixed V19 fit against frozen dev3 only. Acceptance remains per-interface force-vector RMSE ≤0.05 eV/Å and separation-force RMSE ≤0.10 eV/Å. A passing fit still requires a separate data-role audit before formal training integration.

## Status

Geometry screening and the read-only environment preflight are complete. DFT labels are pending; result values will be added only after convergence and archive verification.
