# Ti3SiC2 48-atom GPAW segmented SCF runner

Status: **V4 PATCH SUBMITTED / V4 OFFLINE TESTS PENDING H / REAL GPAW RESTART UNVALIDATED / NO DFT STARTED**

This package handles exactly one fixed-geometry Ti3SiC2 single point.
No historical inputs, results or mother scripts were modified.
Branch base: H revision 3bf6a926ecd0b928c1e3b444f769468d0ee25fe4.

## Scientific provenance (locked)

- Final C science adjudication: f84b163f7346b5d35e4c01cfcfed4a4248f6c8a6.
- H correction: 3bf6a926ecd0b928c1e3b444f769468d0ee25fe4.
- Input: NEW_DFT_GEOMETRY_CANDIDATES/structures/TSC48_static_perturb_std0p025_seed2026101603_UNLABELED.extxyz
- Input SHA-256: da3820c477023e0bb423f87550ed82ae0c2a571f53d294c055f78d7d052d5e8b
- Approved H config SHA-256: ca6f365e8b0b8ec847e329d0312ebb51381b21e996050949494f2d3e0669dfe7
- V3 runner SHA-256: compute from this v3 branch and record independently; old V1/V2 hashes are no longer valid.
- Fixed 48 atoms Ti24Si8C16, fixed cell, PBE/PW500, integer standard Monkhorst-Pack (10,10,4)
  with no explicit Gamma override. Fermi 0.10 eV, spin unpolarized,
  PPCG niter=2, Mixer(.05,8,100), energy/density/eigenstates 1e-5,
  maxiter 160 and four MPI ranks.

## V2 correction after H audit (commit 25eb20bbd064ba5a410bb1bd95a6588c5f224366)

- The V1 MPI stop decision was rank-local. V2 uses shared_decision:
  only rank 0 evaluates wallclock, STOP_REQUESTED, iteration limit and free disk;
  one NumPy int32 vector is broadcast to all ranks before collective GPW writing.
  collective_stage_clear also broadcasts the staged-file status.
- The append-only segment_ledger.jsonl records SEGMENT_STARTED, SCF_SAMPLE,
  CHECKPOINT_COMMITTED, SEGMENT_PAUSED, SEGMENT_CONVERGED with timestamps,
  exact source checkpoint SHA, checkpoint SHA, SCF log path, cumulative
  elapsed time, and cumulative *logged* SCF iterations. These are accounting
  quantities, NOT proof that GPAW resumes the same SCF iteration/mixer.
- Slot manifests include cumulative counters and source checkpoint identity.
  Per-segment GPAW logs preserve the SCF trend; optional numeric SCF context
  fields are emitted only if observed (no guessed convergence numbers).
- ACTIVE_RUN.json prevents concurrent process writes. After a crash it
  intentionally remains: H must verify the previous job is gone before
  any manual recovery action. Never blindly delete an active-job lock.
- On Python-level collective writer exceptions the code attempts world.abort(9).
  This does not guarantee recovery from every MPI/GPW I/O failure, and
  does not replace an authorized real four-rank fault-injection test.

IMPORTANT STORAGE ACCOUNTING: two retained mode=all GPW slots plus a third
temporary staging file may simultaneously exist during rotation. If each
is roughly 11.50 GB, transient GPW bytes alone could exceed 34 GB BEFORE
scratch; prior 23 GB estimate covered only two files. No actual write-size
or scratch peak has yet been measured. This is a safety observation, not a
new arbitrary disk threshold.

V2 is not production validated: H must re-run syntax and FIVE offline
tests before any execution, followed by independently authorized real
GPAW checkpoint/restart testing and genuine resource evidence.

## V3 correction after H V2 acceptance (report commit 115069e0e015395badbd2d77fa946ad6d1df127b)

V2 had two review findings:
1. the cumulative elapsed time stored in a checkpoint manifest was sampled before
   checkpoint writing, so the next segment inherited an undercount;
2. the SCF JSON sampler probed fields that are not exposed by GPAW 26.7.0
   `SCFContext`.

V3 fixes only these two items:
- `save_checkpoint` now measures elapsed time from the segment start through
  GPW publication on rank 0, stores `cumulative_elapsed_seconds` and
  `checkpoint_write_elapsed_seconds` in the slot manifest/status/ledger, and
  broadcasts the post-publication cumulative value to all ranks. The retained
  fallback slot is not overwritten until the new GPW has been successfully
  published.
- `scf_context_sample` now uses the GPAW 26.7.0 `SCFContext` interface:
  `ctx.niter`, `ctx.ham.e_total_extrapolated`, `ctx.dens.error`,
  `ctx.wfs.eigensolver.error`, `ctx.eig_error`, and
  `ctx.wfs.nvalence`. The exported values are explicitly converted to the
  same units used by GPAW's convergence criteria; the authoritative detailed
  trend remains the GPAW segment log.
- Offline tests increase from 5 to 7: one verifies the GPAW 26.7 context mapping
  and one uses a fake checkpoint writer with a deliberate delay to verify that
  cumulative accounting includes checkpoint write/publication time. These are
  not real DFT or MPI-write tests.

GPAW's documented 26.7 SCF loop yields the context to the `scf_step` hook
before convergence checking, and its context exposes the fields used above.
The built-in criteria convert energy/eigenstate/density quantities using the
same `Ha` and valence-electron normalization. See the upstream GPAW source
and convergence-criteria documentation. Real execution still requires H
validation in the actual installed 26.7.0 environment.

V3 does NOT change the scientific input, PBE/PW500, 10x10x4 k-grid, smearing,
solver, mixer, convergence thresholds, atom coordinates, or MPI rank count.

V3 remains unvalidated for actual four-rank GPW write/read, interruption,
persistent storage, scratch peak, and wallclock. No DFT has been started.

## V4 fixes after H V3 static acceptance failed

H V3 reported CHECKPOINT_TIMING_ACCOUNTING=FAIL and
SCF_CONTEXT_MAPPING=FAIL, despite seven offline tests passing. V4 fixes the
specific implementation/test mismatches, without modifying physics:

1. **GPW write-only timing**: record rank-0 write_finished immediately
   when calc.write(mode='all') returns, BEFORE the subsequent MPI barrier.
   Only this interval is checkpoint_write_elapsed_seconds. Barrier,
   SHA-256, fsync and manifest publication are not included.

2. **Cumulative elapsed / safe publication**: fsync + SHA the staged GPW,
   atomically publish it, write an UNFINALIZED slot manifest, append the
   CHECKPOINT_PUBLISHING ledger record, and write the interim status.
   ONLY AFTER all those operations, sample cumulative elapsed and atomically
   rewrite the manifest as checkpoint_finalized=true. verified_slots()
   rejects absent, incomplete, unfinalized or mismatched metadata.
   The final manifest fsync and later MPI broadcast take additional time;
   save_checkpoint() returns a later post-publication elapsed snapshot.
   A manifest cannot encode its own future fsync duration. A graceful
   segment stop records a later end-of-segment cumulative timestamp.

3. **SCF context counter**: record actual ctx.niter as
   scf_iteration_from_context, separate from local_iteration.
   The unit test intentionally makes these two numbers different.

4. **Seven regression tests**: a deliberately slow mock MPI barrier proves
   it is excluded from GPW-only timing; cumulative time includes GPW
   publication and metadata operations; unfinalized GPWs are excluded
   from restart even if SHA and size match. These tests use mock I/O.

H must repeat py_compile, all seven offline tests, exact SHA comparison,
and GPAW 26.7.0 compatibility review. No GPAW scientific computation
has run. RESOURCE_READY=NO, DFT_EXECUTION_READY=NO.

## H offline acceptance (NO GPAW needed)

From a fresh checkout of branch chatgpt-tsc48-segmented-runner-v4-20261016:

~~~bash
ROOT=research/mace-v12-transfer/coordination/reports/chatgpt-nep-ace-dft-repair-20261016
RUNNER="$ROOT/TSC48_SEGMENTED_RUNNER_20261016/tsc48_segmented.py"
INPUT="$ROOT/NEW_DFT_GEOMETRY_CANDIDATES/structures/TSC48_static_perturb_std0p025_seed2026101603_UNLABELED.extxyz"
HCFG="$ROOT/TSC48_PREFLIGHT_20261016/H_REVISION_C_25bb7406/calculation_config_preflight.json"
sha256sum "$RUNNER" "$INPUT" "$HCFG"
python -m py_compile "$RUNNER"
python -m unittest discover -s "$ROOT/TSC48_SEGMENTED_RUNNER_20261016" -p 'test_offline.py' -v
python "$RUNNER" verify --input "$INPUT" --h-config "$HCFG" \
    --workdir /ACTUAL_VERIFIED_PERSISTENT_ROOT/TSC48_singlepoint
~~~

Verify checks the source and config SHAs and does NOT load GPAW.
If the specified workdir does not exist, it checks the input/config only.
The seven included unittest cases use fake checkpoint bytes and a fake MPI broadcast, NOT DFT or real GPAW MPI.

## Minimum real H acceptance before --execute

1. Confirm the actual GPAW 26.7.0 / ASE 3.29.0 / 4 MPI rank environment.
   Review the code and callback integration with the live GPAW version.
2. **Actually validate** mode="all" wavefunction checkpoint writing on four MPI ranks,
   its SHA integrity, dual-slot fallback, real GPAW readback and SCF continuation,
   plus the response to an incomplete .writing.gpw file.  These operations are
   NOT validated by the fake-file tests.  Any extra scientific DFT used in
   such a test requires its own scope approval.
3. Confirm a persistent filesystem across sessions and independent restore
   route. Prove scratch/checkpoint worst-case disk budget, measured write
   reserve, and externally enforced wallclock allowance.
4. Copy resource_clearance_TEMPLATE.json to a site-specific evidence record.
   Only mark true after measuring/validating each item.  The code verifies
   the assertions but cannot independently attest to cloud persistence.
   The shipped template intentionally has resource_ready=false.
5. H must approve a single active process per workdir.  Never start/resume
   concurrently.  Do not change the pinned program after a checkpoint has
   been made, as its runner SHA is part of the immutable lock.

## Segmented interface (EXAMPLES, not a command to run now)

Requires a verified persistent directory and a separately issued, honest
resource-clearance document. Complete all above criteria before executing.

~~~bash
source /workspace/.setup/runtime.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
JOB=/ACTUAL_VERIFIED_PERSISTENT_ROOT/TSC48_singlepoint
CLEAR=/ACTUAL_PATH/resource_clearance_APPROVED.json

mpiexec -n 4 /workspace/.venvs/gpaw-mpi/bin/python "$RUNNER" start \
  --input "$INPUT" --h-config "$HCFG" --workdir "$JOB" \
  --clearance "$CLEAR" --segment-iters 20 --checkpoint-every 20 \
  --segment-seconds 3600 --margin-seconds 600 --execute

# Only after first segment finished and no previous process remains:
mpiexec -n 4 /workspace/.venvs/gpaw-mpi/bin/python "$RUNNER" resume \
  --input "$INPUT" --h-config "$HCFG" --workdir "$JOB" \
  --clearance "$CLEAR" --segment-iters 20 --checkpoint-every 20 \
  --segment-seconds 3600 --margin-seconds 600 --execute

# From another terminal request a safe pause at the next SCF callback:
python "$RUNNER" stop --input "$INPUT" --h-config "$HCFG" --workdir "$JOB"

# Independent offline manifest/byte verification:
python "$RUNNER" verify --input "$INPUT" --h-config "$HCFG" --workdir "$JOB"
~~~

Important: stop is a **request**; it is not instant kill. The two slots
A/B rotate, leaving the last verified checkpoint on the other slot.
A staged .writing.gpw is never a candidate for restart. Do not manually
delete old slots or stale files without explicit H incident review.

SCF convergence is required for RESULT.json. That output remains
DATA_ROLE=HOLD and cannot be registered as a validated training label
merely because SCF converged.  The disturbed candidate's own k-mesh
convergence has not been measured. No subsequent k mesh is authorized.

## Honest current status

RUNNER_V4_PATCH_SUBMITTED=YES
PYTHON_COMPILE_PASS=AWAITING_H_V4
OFFLINE_TESTS_7_PASS=AWAITING_H_V4
CHECKPOINT_WRITE_READY=UNVERIFIED
REAL_GPAW_RESTART_VALIDATED=NO
PERSISTENT_STORAGE_READY=NO
RESOURCE_READY=NO
DFT_EXECUTION_READY=NO
DFT_STARTED=NO

The code itself does not create a persistent volume, certify wallclock
limits, perform a physical DFT restart test, or launch an unattended job.
