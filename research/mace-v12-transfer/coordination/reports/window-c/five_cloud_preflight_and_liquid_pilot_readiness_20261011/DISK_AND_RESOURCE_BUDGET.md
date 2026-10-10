# C read-only disk/runtime review (2026-10-11, Asia/Shanghai)

This is preparation only. No MPI job, DFT, MD, installation, filesystem cleanup, or A/B runtime/checkpoint access was performed. Published compact result logs and report metadata below are allowed public numerical-control material; their roles stay unchanged. Existing 30 GB guidance is retained. D–H machine measurements cannot be inferred from C's machine.

## Public evidence

At original handoff `bb8ad7e4b725b634de2e2198d141df2bfe9d78b9`:

- `research/mace-v12-transfer/coordination/reports/window-c/five_cloud_dft_parallel_plan_20261010/RESOURCE_AND_TIME_BUDGET.md` gives >=30 GB working-disk guidance and >=16 GiB available memory for the 64-atom pilot. These are screening suggestions, not measured size requirements.
- `research/mace-v12-transfer/coordination/reports/window-a/v19_Ti3SiC2_bulk_mesh_comparison/local_checkpoint_cleanup.json` records a completed **48-atom** Ti3SiC2 k8x8x4 restart at **7,361,378,540 bytes** (7.361 GB / 6.856 GiB), historical SHA `3a563ea6def61886af7d6be7547b591fad510e6778daff5b75c678f61d0563dc`, compact archive PASS. This is a published historical metadata record, not an inspection of or instruction to remove any current checkpoint.
- Its public archived `periodic_interface_v4/pbe_interface_v19_pure_phase_controls/calculations/Ti3SiC2_baseline_k8x8x4_v19/gpaw.log` records 128 IBZ points, 14,526 plane-wave coefficients, 235 bands, complex128, wavefunctions 6,977,597,440 bytes and projectors 2,167,508,992 bytes. The archived log is 18,463 bytes; six compact files total about 27 KB. Therefore compact archive sizes do not bound restart sizes.
- Public `periodic_interface_v4/pbe_interface_v19_pure_phase_controls/calculations/Ag_baseline_k8x8x8_v19/gpaw.log` records **32 Ag**, 544 valence electrons (17/Ag), 331 bands, 20 IBZ points, 13,930 PW coefficients, complex128, wavefunctions 1,468,898,560 bytes and projectors 208,574,720 bytes. The log is 14,625 bytes; compact files total about 21 KB. No published checkpoint byte count for this Ag job was identified.
- `coordination/reports/window-b/v19_pure_phase_runner_and_bulk_mesh_preparation/run_pw_reference.py` historically writes `mode="all"` every 20 SCF iterations, checks a legacy byte reserve, and keeps state.gpw local. This is **not** a validated disk policy for the new liquid pilot, and its scientific/ownership logic must not be transplanted or run here.
- `coordination/install_parallel_dft_environment.sh` installation path includes a H2 PW-PBE smoke DFT. Do not execute this installation path in the current no-science task; already installed software is reused. A historical `--check` check is narrower but also assumes a particular /workspace layout and does not verify D–H persistence or all runtime variables.

No sealed structures/labels or live A/B state files were read. At current `origin/main` 684e9a887547404d959e2a7b679c9596ea134a8a there are no `reports/window-d` through `reports/window-h` directories; raw machine preflight evidence may need owner publication. This absence is not proof that the hosts are missing software.

## Decimal GB, GiB, quota and free bytes

- 30 GB = 30,000,000,000 bytes = 27.940 GiB.
- 28.8–29.0 GB = 28,800,000,000–29,000,000,000 bytes = 26.822–27.008 GiB, if the original report used decimal GB.
- If owners actually reported GiB or `df -h` values, the comparison changes. Capture byte count with `statvfs(DFT_WORK_ROOT).f_bavail * f_frsize`, filesystem type/device/mount, filesystem/inode/quota limits and observation time. Effective writable budget is min(f_bavail bytes, remaining quota). Do not measure only the repo filesystem when work root is elsewhere.
- A single host's 28.8–29.0 GB is not five pooled disks. Root path labels alone do not establish separate physical storage or reserved CPU resources.

Current conclusion: **DISK_PILOT_FIT = UNVERIFIED**, and **30 GB guidance remains**. All five owner feedbacks are now received: D/E/F/G/H free bytes are 28797075456 / 28800167936 / 28874539008 / 29018927104 / 28896153600. These are OWNER_REPORTED snapshots, not archived reports/C remote checks. All are below 30 decimal GB. Increasing persistent capacity to >=30 GB resolves only that screening gap; denser meshes and retention can still need more.

## Explicit estimate, not validated sufficiency

Use frozen cell/PAW/method to calculate actual `Nbands`, `Nspin`, `NIBZ`, PW dimensions, FFT dimensions, projector arrays and output settings. A first-order wavefunction component is:

`W_bytes = 16 * Nspin * Nbands * NIBZ * Npw` (complex128).

`Npw ≈ V/(6*pi^2) * (Ecut/3.80998212 eV·Å²)^(3/2)` is a continuum count, not the exact GPAW count. Gamma real/complex storage and time reversal must be verified, not assumed.

The following **engineering illustration only** uses rho=9.3 g/cm3, **which is NOT a verified or approved 1250 K density and MUST NOT set any coordinate box**. At 64 Ag and PW500, it implies V≈1232.65 Å³, Npw≈31,294. Doubling the historical 32Ag 331-band count gives an illustrative 662 bands; actual new bands are UNVERIFIED. This assumes one spin channel and complex128. With raw worst-case IBZ counts 1/8/27:

| mesh scenario | wavefunction component GB | illustrative checkpoint C GB =1.25W |
|---|---:|---:|
| Gamma, 1 point | 0.331 | 0.414 |
| 2x2x2, raw 8 points | 2.652 | 3.315 |
| 3x3x3, raw 27 points | 8.950 | 11.187 |
| 3x3x3, if independently verified time reversal gives 14 points | 4.641 | 5.801 |

The 1.25 overhead multiplier is an unvalidated budgeting hypothesis, not a measured GPAW constant. It can miss densities, projectors, duplicated buffers, rank-local files, larger band counts and writers. GPAW RAM also exceeds saved wavefunction size; total memory and per-rank RSS remain UNVERIFIED.

Budget the **current + previous + atomic temporary** checkpoint (three concurrent copies), all retained checkpoints for earlier meshes/tasks, inputs/outputs/logs and filesystem reserve:

`B = retained_checkpoint_bytes + 3*Cmax + 1,000,000,000 miscellaneous bytes`

`required_free_bytes = B + max(5,000,000,000 bytes, 0.20*B)`.

The 1 GB miscellaneous allowance and reserve are proposed risk controls requiring approval, not modifications to earlier host rules. One active MPI job per machine is assumed. Do not delete completed checkpoint copies to force a fit.

| illustrative retention state | required free GB |
|---|---:|
| two retained Gamma copies + active 2x2x2 checkpoint triple | 16.773 |
| two Gamma + two completed 2x2x2 copies + active second 2x2x2 triple | 23.402 |
| same retained Gamma/2x2x2 files + active 3x3x3 triple (14 IBZ verified) | 31.032 |
| same retained files + active 3x3x3 triple (27 raw points) | 50.422 |

Thus a **limited Gamma/2x2x2-only exception might be feasible** at 28.8–29.0 GB after actual bytes/configuration/retention evidence and a user-approved disk exception. It is not currently demonstrated. A 3x3x3 escalation can fail even if the ordinary 30 GB suggestion is met. Escalation should be a separate capacity/budget/authorization gate; never silently purge checkpoints or change the scientific mesh to fit storage. All six final pilot labels and D's numerical-control checkpoints must enter the cumulative per-host retention ledger. D's 2 pilot structures are not the same structure as E/F's 4 structures; only the designated numerical comparison shares one exact frozen geometry hash.

## Durable independent work root (plan only)

For each X in d/e/f/g/h, propose `DFT_WORK_ROOT=/workspace/dft-work/cloud-X` **only if** the cloud provider/owner confirms that this filesystem is persistent across normal session restart, has adequate quota, is not an A/B or other cloud's mounted run directory, and has separate write ownership. If /workspace is ephemeral, choose a verified persistent attached volume before authorization. `/tmp` is not a durable choice. No root was created or configured by C.

Keep config outside run and source repository, for example `/workspace/dft-config/cloud-X/runtime.sh`; use private branch `cloud-X-b1-dft`. Runs use `${DFT_WORK_ROOT}/<task_id>/<numerical_run_id>/{input,checkpoint,logs,output}`. The numerical_run_id freezes exact method/kmesh and avoids sharing a checkpoint across mesh variants. Approved input packages can be mirrored read-only by SHA; checkpoint files and rank scratch never shared across hosts. Record realpath, mount source/type, filesystem identity hash, quota, owner, provider's independent instance/allocation attestation, and persistence statement. Different directory names on the same/shared volume are not enough to prove host-resource independence.

For isolation checks, rely on owner manifests/provider attestations and allowed mount/allocation metadata; do not inspect or modify A/B working directories/processes. Compare only published path namespaces/manifests. Mark any underlying physical separation without evidence UNVERIFIED.

## F/G/H environment runtime checklist

Versions being installed does not prove fresh-shell runtime works. For each owner:

1. Publish the exact interpreter and launcher path/version, verified environment root, `runtime.sh` existence/path/SHA, resolved `GPAW_SETUP_PATH`, and Ag/Ti/Si/C PBE setup file SHA matched to the project manifest. Keep private home paths/credentials out of reports.
2. Inspect runtime.sh as text before sourcing. It should only configure the validated environment; no package installer, scientific launcher, background job, cleanup, or automatic checkpoint restart. Missing file is CONFIGURATION_UNVERIFIED, not INSTALL_REQUIRED.
3. Configure the existing GPAW environment PATH and library locations, `GPAW_MPI_BACKEND=cgpaw` if this is the audited build, `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`. Resolve `GPAW_SETUP_PATH` from the installed gpaw_data.datapath() **before importing GPAW** and fail if directory/setup SHA differs. Do not append an unverified second setup directory that changes lookup precedence.
4. Check `gpaw info` for MPI/ScaLAPACK build and linked libraries. Compare OpenMPI launcher/runtime/ABI to `_gpaw` linkage and provider environment. Merely having mpiexec 5.0.7 is insufficient. Do not recompile/reinstall working components.
5. Under any separately approved non-scientific runtime sanity check, a short `mpiexec --bind-to core --map-by core -n 4 <verified-python> -c 'from gpaw.mpi import world; assert world.size == 4'` verifies rank initialization only; it must not construct Atoms/calculator or request energy/forces. This review did not run it. Do not add oversubscribe or root bypass flags automatically.
6. Recheck effective CPU allocation: sched_getaffinity/cpuset count, CPU quota (`cpu.max` or v1 quota/period), actual physical topology and provider reserved allocation. Effective schedulable capacity is bounded by affinity, cpuset and quota, not /proc's visible logical CPU count. Check memory limit (`memory.max - memory.current`) as well as MemAvailable; do not use host-wide free RAM if a cgroup quota is lower.

Update from owner messages: G reports runtime.sh missing, with existing /workspace/.setup/activate.sh and previously matching GPAW_SETUP_PATH; H reports /workspace/.setup/runtime.sh, matching setup path, OpenMPI 5.0.7 and linked MPI/ScaLAPACK with no missing libs. These are OWNER_REPORTED, not C remote verification. F reports /workspace/.setup/onboarding/runtime.sh loading, four setup hashes matching, resolved libraries and historical 4-rank PASS; DFT_WORK_ROOT remains unset. E reports /workspace/.setup/e/runtime.sh syntax/setup PASS and work-root BLOCKED. D reports /workspace/onboarding/runtime.sh syntax PASS and the required GPAW_MPI_BACKEND=cgpaw for a successful 4-rank probe. Full script SHA/content, persistence and original reports are still missing. See HOST_READINESS_MATRIX.csv for current exact bytes.

## Proposed runtime/disk stops after explicit authorization

- Before each task and before any checkpoint write, recompute remaining bytes/inodes/quota and predicted additional bytes needed for temporary/current/previous writes plus reserve. Do not start or write if the risk margin would be violated; preserve existing state and report. Implement/write behavior itself must be approved ahead of stage 2; C does not change owner scripts now.
- Never delete/truncate an existing checkpoint or logs to make room. Atomic replacement policy retains the last fully verified recoverable state; write/fsync/rename semantics and hash/checkpoint metadata must be reviewed before implementation. Stop if rename/write fails; do not continue silently without checkpoint coverage.
- Unexpected checkpoint size over budget, a denser kmesh, extra bands/spin, large rank-local files, changed mount/quota/persistence or cross-cloud path collision = stop corresponding task and request re-budget/approval.
- Memory reserve and effective 4-core allocation must be met; one 4-rank job per independent host. No oversubscription or cross-cloud sharing to rescue resource gaps.
- First approved liquid calculation beyond 8 h stops further batch submission for re-estimation; per-run proposed 12 h/160 SCF limit preserves checkpoint and reports. No automatic restart or silent parameter/geometric changes.
- This note recommends separate approval of (a) real geometry-only sampling, (b) frozen-coordinate numerical mesh pilot and checkpoint policy, (c) six final labels, (d) denser escalation/remaining24. None is approved by this preparation report.
