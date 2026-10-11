# F DFT preparation report

Captured: 2026-10-11 (Asia/Shanghai)
Owner: F
Branch at capture: `cloud-f-b1-dft`
Scope: prepare the six F-window DFT work items without generating structures or running calculations.

## Status

The GPAW route is available for future approved work. The QE route is unavailable in this environment. The mounted workspace has free space, but it is a task-temporary overlay and no DFT work root exists; storage is therefore **not accepted for checkpointed DFT work**. All six structures remain `STRUCTURE_NOT_READY`; the two Ag64 liquid tasks have no approved candidate input, source hash, or task-specific k-point record.

| Item | Finding | Status |
|---|---|---|
| Python / GPAW / ASE | Python 3.12.14, GPAW 26.7.0, ASE 3.29.0 import successfully in the existing F runtime | PASS |
| GPAW build | GPAW reports MPI enabled, OpenMP enabled, ScaLAPACK enabled, LibXC 5.2.3; runtime selects `GPAW_MPI_BACKEND=cgpaw` | PASS for import/build inspection; no MPI calculation launched |
| QE | `pw.x` is not installed; no `.UPF`/`.upf` file is present in the repository | NOT AVAILABLE |
| PAW data | gpaw-data 1.2.1; Ag, Ti, Si and C PBE files exist and match the frozen project SHA-256 baselines below | PASS |
| PAW licensing metadata | gpaw-data distribution metadata declares GPLv3-or-later. GPAW's installed distribution metadata has no `License` field; no separate redistribution review was performed | PARTIAL |
| Workspace storage | `/workspace` is an `overlay` filesystem with 28,869,459,968 available bytes and 2,036,353 free inodes at capture; this is 1,130,540,032 bytes below the handoff's 30 GB recommendation | MEASURED, TEMPORARY |
| DFT work root | `/workspace/.work-f/dft` and `/workspace/dft-work/cloud-f/ag-ti-si-c-b1` do not exist | NOT READY |
| Checkpoint budget | No F Ag64 candidate or representative GPAW checkpoint exists to measure; temporary storage is not durable evidence | UNVERIFIED |

PAW file hashes:

| Element | File | SHA-256 |
|---|---|---|
| Ag | `Ag.PBE.gz` | `bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97` |
| Ti | `Ti.PBE.gz` | `bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f` |
| Si | `Si.PBE.gz` | `2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0` |
| C | `C.PBE.gz` | `fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983` |

`GPAW_SETUP_PATH` matches `gpaw_data.datapath()` at `/workspace/.venvs/gpaw-mpi/lib/python3.12/site-packages/gpaw_data/setups`. The existing runtime script is `/workspace/.setup/onboarding/runtime.sh`, SHA-256 `299f595beb8233cc7b22601a01d1da7e63380ce4cc044d6cd67d8bf979a9ed00`. GPAW/ASE checks used that runtime's Python 3.12.14 executable. No software or pseudopotential was installed or replaced.

## Six-item F window checklist

| # | Work item | Prepared result | Remaining gate |
|---:|---|---|---|
| 1 | DFT environment: ASE, GPAW/QE, pseudopotentials, storage | Environment and PAW evidence above; QE absence and temporary-storage limitation recorded | Durable approved work root and measured checkpoint budget; QE only if C selects that backend |
| 2 | Ag single-point flow | `AG64_SINGLE_POINT_FLOW.md` and `AG64_INPUT_MANIFEST.template.json` | C-reviewed liquid candidate, density/source lineage, exact input SHA, k-point pilot record, release |
| 3 | Ti₃SiC₂/Ag interface DFT preparation | `INTERFACE_PREPARATION.md` and `INTERFACE_INPUT_MANIFEST.template.json` | There is no mixed Ag/Ti₃SiC₂ interface task among F's six assigned IDs; C must assign a task and freeze geometry/source/method before use |
| 4 | Input templates | Ag64 and interface JSON manifest templates | Replace placeholders only from reviewed input evidence; no coordinates are included |
| 5 | Energy, force and stress format | `RESULT_FORMAT.template.json` and output conventions in `AG64_SINGLE_POINT_FLOW.md` | Populate only from an approved calculation; none was run |
| 6 | DFT result archive format | `ARCHIVE_LAYOUT.md` | Create a task archive only after a candidate, release and durable storage are available |

The six F proposal IDs in `HANDOFF_F.md` remain unchanged:

- `C-F-AGL-1050K-rho103-seed01`
- `C-F-AGL-1250K-rho100-seed02`
- `C-F-AGL-1250K-rho103-seed01`
- `C-F-AGL-1600K-rho097-seed01`
- `C-F-TSC-PF01-CFG03`
- `C-F-AG-PF01-CFG03`

Every item still lacks frozen coordinates and a runnable input SHA. The two 1250 K liquid tasks are candidates in the C handoff, not structures. The DFT template requires an existing input path and an explicit execution enable flag; this preparation did not invoke it.

## Stop state

No atom structure was generated or modified. No GPAW or QE calculation, MD/AIMD/LAMMPS/TTM simulation, or training run was started. Keep all six F tasks on HOLD until reviewed structures and an explicit release are available. For liquid Ag, use the same frozen 64-atom geometry for the Γ and 2×2×2 comparison; apply the project's pre-registered energy and force-difference criteria before choosing a mesh. Do not use Γ by default.
