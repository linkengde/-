# H minimal two-file relay implementation — offline only

## Outcomes

- Branch: `h-tsc48-minimal-relay-20261016`.
- Actual Git push succeeded using existing platform Git identity. No new token requested.
- Repository API inspection: `linkengde/-` is PUBLIC. No scientific GPW was uploaded.
- 5 MiB random synthetic file (5 x 1 MiB parts) really uploaded in immutable commit
  `bdf0f292edc1d75bbfe41655f70373838ca9a4ab`.
- A fresh Python process fetched the manifest and downloader scripts from that commit,
  streamed five real remote parts, and reconstructed SHA
  `ce41df6ec4ca386de691755fa545f4b2f9255b95b186795c8bf40910cdd79022`.
- Remote readback: 1.818 seconds, process maximum RSS 18,716 KiB. Payload 5,242,880 bytes;
  protocol overhead and upload peak were not instrumented. No full local chunk set at receiver.
- This reader was on H1 in a different process, NOT an independent H2 window.
  `CROSS_WINDOW_DOWNLOAD_PASS=UNVERIFIED`. Execute `H2_SYNTHETIC_READBACK.sh` in H2
  and return its real evidence; this H session has no H2 execution interface.
- Four actual GPAW MPI ranks exercised the adapter with a FAKE calculator writing
  synthetic bytes: two files retained and third write blocked on all ranks. No DFT.
- 14 offline tests pass: corruption, partial write, pointer publication failure,
  two-file cap, receipt refusal and predecessor retirement, plus imported chunk tests.
- Input/config verification and Python compilation pass. No scientific calculation started.

## Implementation and limitations

`frozen_v4.py` is byte-for-byte V4 from commit
`5b67aa3002b1e3c197d7d511e9e075586c942095`, SHA
`23250241f51004380850233cabbbed443be38087c47b286977bdead680ed5c95`.
The original V4 path is unchanged. `relay_runner.py` overrides only checkpoint storage,
selection and default artificial segment/checkpoint iteration limits. Calculator construction,
scientific parameters, PAW/input/config checks, MPI rank-0 decisions and site resource
clearance remain V4. Defaults use very large iteration cutoffs, so normal SCF is continuous;
external wallclock/stop request remains the planned handoff trigger.

On collective write all ranks write one `checkpoint.writing.gpw`, preserving the old valid
GPW. After fsync/hash, rename to a unique generation, publish metadata, then atomically
publish `current.json` as commit pointer. A publication failure preserves the old pointer;
an orphan candidate blocks further writes for incident review. Exactly two GPW files are
allowed; a third write is refused until the previous generation is safely retired.
This is bounded storage, NOT guaranteed DFT warm-restart reliability.

`retire_previous` is explicit maintenance only when no ACTIVE_RUN lock remains. It requires
an independent-reader receipt bound to the current complete GPW SHA, byte count and remote
immutable identifier, verifies the previous file SHA, and never deletes the current GPW.
The receipt is an operator attestation: JSON flags cannot prove a remote reader actually
ran or guarantee retention. No production receipt exists and no scientific files were deleted.
Crash recovery of orphan stages, stale locks and ambiguous generations is manual/fail-closed.
There is no automatic cross-window scheduling or global distributed execution lock.

Minimal handoff: latest full GPW; current.json; lock.json; ordered chunk manifest with all
SHA/size fields; source/sequence/cumulative SCF accounting; pinned V4/adapter/input/config
commits and hashes; GPAW/ASE/MPI/PAW environment evidence. Small audit ledger excerpts may
be archived. Do not carry ACTIVE_RUN or STOP_REQUESTED into a fresh receiver. Wavefunction
warm restart does not preserve the SCF mixer history bit-for-bit. Real four-rank GPW write,
fresh-process readback and genuine continued SCF are still UNVALIDATED.

## Disk budget (not measured scientific usage)

Observed workspace free before test: 28,826,624,000 B (28.83 decimal GB).
After adding synthetic parts and Git objects: 28,815,618,048 B; the roughly 11 MB change
already illustrates worktree plus Git-object storage overhead. `/tmp` is tmpfs: receiver
measurements there are not proof of a persistent disk or DFT RAM safety.
For old/current GPW S_old and new GPW S_new, required local scientific peak is
S_old + S_new + measured scratch + metadata + reserved margin, on top of baseline use.
At the UNMEASURED 11.5 decimal GB/file estimate, two files consume 23 GB, leaving about
5.82 GB. A single 32 MiB transfer chunk consumes 33,554,432 B, not another whole GPW.
If '11.5' means GiB instead, two files consume 24.70 decimal GB and leave only about
4.12 GB. Actual GPW size, write duration, DFT peak RAM, scratch and write reserve are UNKNOWN.
Do not declare the two-file budget safe until measured. Checksum validity alone is not a
proof that wavefunctions can be used for successful warm restart.

The inherited split tool still writes the complete local chunk set; use it ONLY for this
5 MiB fixture here. A production bounded upload layer remains unfinished. Streaming
receiver is implemented. A per-chunk upload layer must account for Git's local object
database as well as worktree chunks; unlinking local files does not reclaim Git objects.
Normal Git history accumulates every multi-GB GPW version even after later deletions.
32 MiB chunks avoid the individual-file limit but do not remove push-size limits, history
bloat or confidentiality problems. No production checkpoint can use this public repository
without explicit publication approval, and ordinary Git is not qualified here for recurring
multi-GB checkpoints. No force push/history rewrite or paid storage was performed.

## Go/no-go

H1_DFT_CAN_START=NO
H2_SCIENTIFIC_CHECKPOINT_AVAILABLE=NO
GIT_ACTUAL_PUSH=PASS
REMOTE_SYNTHETIC_READBACK=PASS
CROSS_WINDOW_DOWNLOAD_PASS=UNVERIFIED
REAL_GPAW_RESTART_VALIDATED=NO
RESOURCE_READY=NO
DFT_EXECUTION_READY=NO

Remaining gates: actual independent H2 synthetic receipt; suitable authorized scientific
checkpoint storage and bounded uploader; measured two-file GPW/scratch/RAM budget; genuine
four-rank mode=all write/read/warm restart evidence; platform external wallclock allowance.
No change to approved science and no request for a 100 GB disk or purchased server.
