# F Codex Cloud temporary-storage substitute acceptance

Date: 2026-10-11 UTC  
Repository branch: `cloud-f-b1-dft`  
Base HEAD before this addendum: `202cae550a6592a75b7431066396a1d88ed448ef`

## Result

- `OWNER = F`
- `WORKSPACE_WRITE_VERIFIED = YES`
- `TASK_LOCAL_SPACE_READY = YES` (28,873,072,640 available bytes; 2,036,714 free inodes at capture)
- `PERSISTENCE_VERIFIED = NO` (Codex Cloud task workspace is temporary by design)
- `GITHUB_UPLOAD_DOWNLOAD_SHA_RESTORE = YES` for the synthetic recovery fixture
- `SAMPLER_ROOT_READY = NO`; `DFT_WORK_ROOT_READY = NO` (both absent and intentionally not created before C's signed path config)
- `SAMPLING_STORAGE_READY = NO` pending C's Codex-specific signed release config and exact roots
- `DFT_STORAGE_READY = NO`; DFT and full sampling continuation remain HOLD
- `WARMUP_RELEASE = HOLD` (`storage_release=false`, `warmup_release=false` in current config SHA `57dc948d45daa3c51f7b3e5491a19e4f407546806531c3afae14a221b054ffc0`)

## Workspace boundary

`/workspace` is a writable overlay with 33,770,192,896 total bytes and 28,873,072,640 bytes available (26.89 GiB). It has 2,036,714 available inodes. The checkout is `/workspace/-` (mode 0700); the managed write roots are `/workspace` and `/tmp`. A new report directory and a 16,415-byte deterministic non-science test file were successfully written under the checkout. The overlay does not provide evidence of persistence after task recreation. The owner sampling and DFT root candidates checked for this report are absent, and no root was created.

The cgroup CPU quota is 4 CPUs (`400000 100000`); memory limit is 34,359,738,368 bytes, with 30,780,076,032 bytes of limit headroom at capture. No GPAW, LAMMPS, or MPI launcher process was present. These are task-time observations, not durable capacity guarantees.

## GitHub branch round-trip

The synthetic fixture was committed and pushed by a normal fast-forward to `cloud-f-b1-dft`. A fresh temporary Git repository then fetched the branch from GitHub at commit `afb1ec415df43dbf6a48212723237d8a999fd00a`, downloaded the fixture from that fetched commit, copied it into a separate restore path, and compared both copies with the source. All three SHA-256 values matched `24f3c1e665f4e0f74ea736d493704b7756871e238d2956761b001aff62921c6b`; byte count was 16,415. This verifies upload, remote fetch, content integrity, and byte restoration for the non-scientific fixture. It does not prove that scientific archive bytes are acceptable under the repository's visibility/data-handling policy or that task-local files persist.

## Ag64 warmup size and time estimate

The sampler writes an initial frame plus a frame every 100 steps, so a 1000-step warmup produces 11 trajectory frames. It writes thermal rows at step 0 and every 20 steps, for 51 rows. The final restart bundle duplicates the current trajectory and thermal file and adds state arrays, state trajectory, JSON state, and integrity metadata. Based on 64-atom coordinate/momentum payloads and these counts, estimate **0.2–1 MiB per task** including the live files and restart snapshots, or **0.4–2 MiB for F's two tasks**. This is a format/count estimate, not a measured run or serialized scientific output. It is comfortably below the 512 MiB proposal cap and GitHub's per-file limit; actual byte counts must be recorded after an authorized warmup.

No Ag64 benchmark is available, and this task did not run the sampler. For planning only, estimate **about 1–5 minutes per task** for 1000 single-process ASE/EAM steps plus startup and checkpoint I/O, with **up to 10 minutes reserved per task**. The currently proposed hard warmup ceiling is 3,600 seconds per task. These timings are unmeasured and must be replaced with observed `wall_s` and `steps_per_s` from the signed warmup. If performance is materially slower, the runner's configured wall stop applies; do not increase it locally.

## Existing software and reports

Existing Python 3.12.14 / GPAW 26.7.0 / ASE 3.29.0 / OpenMPI 5.0.7 / gpaw-data 1.2.1 and PAW SHA evidence remain in the previously submitted F report; this addendum does not reinstall or alter them. `run_sampling.py` SHA is `b78b734798c9c075826e6611e3d2f7f156f0b05ff7c87159812ac76ed0e918a8`. The prior F resource/runtime report remains at `../REPORT.md` with its existing checksums intact.

Supplementary provenance evidence for the attached Ag release notes and paper is recorded in [AG_POTENTIAL_PROVENANCE_SUPPLEMENT.md](AG_POTENTIAL_PROVENANCE_SUPPLEMENT.md). The supplement records source identity only; it does not change the release config or the separate redistribution hold.

## Required next step

C should publish a Codex-specific config and SHA that accepts the temporary workspace plus tested GitHub archive/recovery path, signs exact F task roots and archive rules, and explicitly releases warmup. Until that arrives, do not run the 1000-step warmup. After the two warmups, report measured size, wall time, restart SHA, and resource observations back to C. DFT and full continuation remain HOLD.
