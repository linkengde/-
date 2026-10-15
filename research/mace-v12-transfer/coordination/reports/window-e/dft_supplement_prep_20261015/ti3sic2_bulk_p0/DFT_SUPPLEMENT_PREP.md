# Ti₃SiC₂ P0 DFT input preparation

**Status: input candidates only; no DFT labels were generated.** This is the first six-point backbone response batch selected from the audited gap list: 2% volume ±, 0.5% basal shear ±, and one Si-site z displacement ±0.02 Å.

## Source and geometry integrity

The parent is copied byte-for-byte from fixed commit `499ddf0870235594444abad780582ea4fcc89d61`:

`research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/inputs/Ti3SiC2_baseline_PROPOSAL.extxyz`

Its raw SHA-256 is `157c1dab02a28e6140dd5191187ccd914fe558950e2ac6c136866bcc7805c3e4`. It has 48 atoms (Ti24Si8C16), PBC `[true,true,true]`, cell `[[6.115,0,0],[-3.0575,5.2957453441,0],[0,0,17.6235]] Å`. The ordered symbols and coordinates match the archived converged `Ti3SiC2_baseline_k10x10x4_v19` result exactly (maximum coordinate difference 0 Å). That result and the Ti/C ±0.02 Å points form the same bulk parent family; this new batch must stay with that family in any later data split.

The generated files preserve atom order, `local_new_id` and space-group kind arrays, use fixed coordinates/cell within each requested geometry, and contain no energy, force or stress labels. The Si representative is `local_new_id=7`, `spacegroup_kinds=2`; all eight Si sites share that archived kind. The shear definition is explicitly `x' = x + gamma*y` with `gamma=±0.005` in the basal plane. Each transformed geometry is an input candidate, not a claim that C has approved that strain convention.

## Existing software and resources

No download was needed: the E-private environment is already present at `/workspace/.venvs/gpaw-mpi`. Read-only checks found Python 3.12.14, GPAW 26.7.0, ASE 3.29.0, gpaw-data 1.2.1, GPAW MPI enabled, OpenMPI 5.0.7 and the project Ag/Ti/Si/C PAW SHA-256 values matching the archived baseline. `/workspace/.setup/e/runtime.sh` (SHA-256 `e32d19eef429939bcb3d0aa2d365d063b6dca16b4b4052b02da226fd7247e43b`) sets `GPAW_SETUP_PATH` from `gpaw_data.datapath()`, `GPAW_MPI_BACKEND=cgpaw`, and OMP/BLAS threads to one. Its E work-root check requires `/mnt/e-work`; that mount is absent. The current Codex workspace has a 4-core quota, 32 GiB memory limit and 22.43 GB free, but no dedicated DFT work root or persistence proof.

## Execution gate

`role_freeze.json` remains unchanged with TRAIN/DEV/TEST empty. The six candidates share the audited 48-atom parent and may not be split across roles. The parameter proposal and candidate manifest are explicitly not frozen, and no new C DFT execution authorization is present. Therefore this package does not start or submit calculations. Before execution, C must freeze the role of the entire parent family and exact inputs; confirm k-point centering, MPI count, wall-time/checkpoint rules and DFT work root; then issue authorization bound to the manifest/structure SHA values.

The broader audit still identifies additional requirements beyond these six points: an independent clean surface family, mixed-interface registry/gap coverage, DFT-labelled liquid Ag and liquid-solid contact/coexistence structures, plus a parent-family held-out TEST set. The six backbone points do not by themselves make NEP/ACE training ready.
