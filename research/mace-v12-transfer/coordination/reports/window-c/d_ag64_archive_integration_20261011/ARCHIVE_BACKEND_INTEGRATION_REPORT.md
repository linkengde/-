# Stage D archive backend integration and bounded execution release

Status: `ARCHIVE_BACKEND_PRODUCTION_PASS=YES` for the GitHub upload, independent fresh-clone verification, and restore interface. A C-signed execution release is present for the single D task only. No MD, force/energy evaluation, DFT, or model training was run in this stage.

## Source and branch identity

- D source checkpoint commit: `145ee8050cf95a3a1988e868c6b935e10538f723`.
- Step-1000 archive SHA-256: `b62a818e58a2cf423e4cb1b9ed0468ae73ad578a75c4b2cb9fb3f28baf75e7b1`.
- Step-1000 outer manifest SHA-256: `ab9f0dffa8ef145eac7792dfe1e14e07c1354df1c2152dafec97bad1dc3a57b1`.
- D archive test commit/head at signing: `c1c9f866426e572e82ab88157654a015802707f8`; its supplied test archive SHA `1f9f5dc0536b2736cf8b05583631bf7c3dcd3f595f73db38a3534d97ff33c924` was independently checked from a fresh D-branch clone, with all listed member hashes passing.
- D branch head was rechecked as the same `c1c9...` immediately before release generation. Any subsequent D-branch update makes this release fail closed.

## Backend and recovery test

`github_archive_backend.py` implements `upload()`, `verify()`, and `restore()`. Upload creates a new hash-addressed commit and refuses to overwrite an existing archive path. It pins the expected branch head, uses a normal fast-forward push, and rejects a moved branch. Verification clones the remote branch into a new directory and checks the archive SHA, member set, manifest, and each internal SHA-256. Restore uses another fresh clone, extracts regular files only into a new path, and verifies every restored member again. The 50 MiB archive/member ceiling is enforced.

The end-to-end recovery test copied the archived D step-1000 checkpoint into a temporary directory, restored its coordinates, velocities, and PCG64 state, then produced a deterministic state-only synthetic step-2000 fixture. It did not create a calculator, evaluate an EAM energy or force, or call an MD integrator. The synthetic archive (`2791a1e8ab1d670481961ad7ceb35d73a31885e584dc484244220c2ae1261361`) was uploaded to the C-owned test path, read back from a fresh GitHub clone, restored, and compared exactly:

- Upload: `PASS` (remote commit `6bc38bdcb13f9ad6cbc38b6c6557878a11b63cd5`).
- Fresh-clone archive and internal hashes: `PASS` (6 payload members).
- Restore: `PASS`.
- Positions: exact match.
- Velocities: exact match.
- PCG64 RNG streams: exact match.
- Original step-1000 archive bytes: unchanged.

The C test archive lives under `backend_roundtrip_tests/` and is explicitly synthetic, not a scientific segment result. D's branch was not modified by this C-side round-trip test.

## Frozen code and execution boundary

The original ASE segment core is copied byte-for-byte. Its hard integrator-owned stop remains global step 2000; `run()` and `irun()` remain disabled. The new entry adapter binds the core, backend, release configuration, source checkpoint, D branch head, and potential hashes before constructing the EAM calculator. Each 100-step checkpoint is uploaded, freshly cloned, verified, restored, and checked before the next integrator step. Any transport or verification failure stops the segment and preserves the original checkpoint and local outputs.

The C-signed execution release is `SIGNED_SEGMENT_EXECUTION_RELEASE.json`; its exact hash is in `SIGNED_SEGMENT_EXECUTION_RELEASE.sha256`. It authorizes only `C-D-AGL-1250K-rho097-seed01`, from step 1000 through step 2000 inclusive, with at most 1000 added steps. It prohibits step 2001, 120000-step continuation, potential/checkpoint/temperature changes, and RNG or velocity reinitialization. The D runtime must still verify the live potential, source archive, unused output root, and exact D branch head before it constructs the calculator.

No D process was started or notified by C. The authorization is available on the C branch for D to retrieve and use with the hash-pinned entry command. It is not authorization for the remaining five warmups, DFT, or training.

## Release hashes

- `CONFIG_SHA256`: `026be260ca65ab55b60c96a407c6181649e5f245bcf2d138993e246d2aa56d9b`
- Production entry adapter: `e855ad26cf215281e04b612163cb9db3c075ca8a6366686cafcec5becde264df`.
- GitHub archive backend: `8fd001a9b1706fc8832e90a77c22783a22a77eab06f50a66bd1d1eeb97c00346`.
- Frozen ASE segment core: `e0bba3702c54bbdcb39f5f7e5641e1e1969a2ea5af8faceee8a0d50316f440dd`.

The configuration permits only the signed step 1000→2000 D segment. It does not start or notify D; it is available for D to fetch from the C branch. D must pass the signed live checks before the adapter constructs the calculator.
