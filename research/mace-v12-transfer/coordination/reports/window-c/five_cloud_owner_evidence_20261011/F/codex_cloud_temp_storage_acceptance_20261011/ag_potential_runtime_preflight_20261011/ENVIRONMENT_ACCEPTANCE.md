# F Ag potential and Codex runtime acceptance

Captured: 2026-10-11 UTC  
Owner: F  
Current repository branch: `cloud-f-b1-dft`  
Previous F branch tip before this report: `4addb5c3ced736589737a005d5c578331291c8d4`

## Required statuses

- `LOCAL_POTENTIAL_PASS = YES`
- `ASE_EAM_PARSE_PASS = YES` (ASE 3.29.0; parser is metadata-only)
- `NODE_PREFLIGHT_PASS = HOLD` (all live node/resource/potential checks pass; C's exact config remains HOLD)
- `REMOTE_ARCHIVE_PASS = YES` (previous F synthetic fixture roundtrip retained, and this report is pushed and re-fetched for SHA verification)
- `BLOCKERS = C config has no scientific release; F host release is false and neither F task is in warmup_task_ids`

## Uploaded ZIP and extracted local file

Input ZIP path: `/workspace/attachments/eaafd8f9-bfa5-44d7-b404-066e2c33a148/WMH2006_Ag_四个非PDF官方文件.zip`  
ZIP size: 395,848 bytes  
ZIP SHA-256: `fd3670c0cf714c60b706c6b2d485ef6134b58b9fdb75aff86ca739ecf846d51c`  
The original ZIP was extracted as supplied, without re-downloading or replacing any data. Its ZIP integrity test passed and all four source-data member hashes match its embedded `SHA256SUMS.txt`.

Extracted potential path: `/workspace/wmh2006-ag-F.RdBqAx/Ag.eam.alloy`  
File size: 780,452 bytes  
SHA-256: `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`  
This matches C's frozen Ag parameter SHA. The potential file itself is local to this task workspace and is excluded from this GitHub report.

## ASE parser result

C's `inspect_eam.py` was taken from the checksum-verified C release branch `cloud-c-codex-ephemeral-warmup-20261011`, commit `f174cf6adbe35101f993a613f725586e27de8bc1`. Its SHA-256 is `d3a437582dcd5bad7e361c6c81f64ac1191161170105b43133ffdd4065526a77`, matching the existing C execution pack. The report identifies ASE `3.29.0`, format `alloy`, element `Ag` / atomic number 47, 10,000 radial and density grid points, and cutoff `5.995011000293092 Å`. It completed the intended metadata parse without evaluating energy or forces, constructing atoms, or running MD.

The parser labels the radial tail `C_TAIL_REVIEW_REQUIRED` because the published setfl cutoff is one grid interval beyond its last tabulated radius. The measured gap is one `dr` (`0.0005995011 Å`), with tail density and pair values near zero. The existing C `SOURCE_RELEASE.json` and sampler review explicitly approve this exact one-bin published setfl tail policy; no table was altered. The raw parser JSON is attached below.

## Current runtime and node preflight

Runtime was sourced from `/workspace/.setup/onboarding/runtime.sh` (SHA-256 `299f595beb8233cc7b22601a01d1da7e63380ce4cc044d6cd67d8bf979a9ed00`). It reports Python 3.12.14, GPAW 26.7.0, ASE 3.29.0, gpaw-data 1.2.1, NumPy 2.5.3, SciPy 1.18.1, `GPAW_MPI_BACKEND=cgpaw`, and OpenMPI 5.0.7. The existing Ag/Ti/Si/C PAW hashes all match the project baselines in the raw helper report.

C's complete Codex release package was fetched read-only from the branch above and all 19 entries in its `SHA256SUMS.txt` passed. C's `node_preflight_codex.py` SHA is recorded in `C_PACKAGE_IDENTITY.json`. It was run separately for both F task IDs, with the exact C config SHA and the extracted local potential path. Both raw reports show:

- CPU quota 4 cores; cgroup memory limit 34,359,738,368 bytes; process RSS about 12 MiB.
- Current `/workspace` free bytes about 28.87 GB and free inodes about 2.04 million; preflight thresholds (2 GiB plus 32 MiB reserve, 10,000 inodes) pass.
- F session fingerprint matches C's recorded F fingerprint.
- The configured roots are absolute, outside the checkout, non-overlapping, on the same filesystem, and the sampler root is absent. `persistence_verified=false` is retained as required for Codex temporary storage.
- The exact local Ag potential SHA check and F's synthetic remote roundtrip evidence pass.
- The only preflight reason is `C_CONFIG_OR_OWNER_RELEASE_HOLD`; the script exits 2 and reports `PREFLIGHT_STATUS=HOLD` for both tasks.

C's config SHA is `207451d2d43279b9a1ac2d661bbca56b60641890781090e9ea14ecb55ebbd34a`. Its attestation says HOLD; `config_status=C_HOLD_NOT_SIGNED_FOR_SCIENTIFIC_EXECUTION`, `release_decision=HOLD`, `storage_release=false`, `warmup_release=false`, and `F_host_release=false`. Its only `warmup_task_ids` entry is D's first task; neither F task is released.

## Remote archive evidence

F's previously qualified synthetic fixture test remains in the parent directory's `../RECOVERY_TEST_RESULTS.md`: a 16,415-byte non-scientific fixture was pushed to `cloud-f-b1-dft`, fetched by a fresh Git repository, restored, and SHA-256 checked. The fixture and report files are metadata/test content only. This report's push and fresh-fetch SHA verification are recorded in the handoff commit.

The GitHub repository is public. This commit contains only hashes, paths, runtime/preflight metadata, and parser output. It does not contain `Ag.eam.alloy`, the PLT files, PDFs, or any other potential payload.

## Stop status

No software was installed. No atoms or structures were generated. No MD, DFT, AIMD, LAMMPS, TTM, or training was run. The current C config is an audited HOLD package, not an authorization to begin warmup. Wait for C's new signed config and exact SHA before starting any science.
