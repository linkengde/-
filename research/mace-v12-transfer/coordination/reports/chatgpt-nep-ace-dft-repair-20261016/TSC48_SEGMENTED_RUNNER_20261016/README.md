# Ti3SiC2 48-atom GPAW segmented SCF runner

Status: **CODE DELIVERED / OFFLINE-TESTED / GPAW RESTART UNVALIDATED / NO DFT STARTED**

This package handles exactly one fixed-geometry Ti3SiC2 single point.
No historical inputs, results or mother scripts were modified.
Branch base: H revision 3bf6a926ecd0b928c1e3b444f769468d0ee25fe4.

## Scientific provenance (locked)

- Final C science adjudication: f84b163f7346b5d35e4c01cfcfed4a4248f6c8a6.
- H correction: 3bf6a926ecd0b928c1e3b444f769468d0ee25fe4.
- Input: NEW_DFT_GEOMETRY_CANDIDATES/structures/TSC48_static_perturb_std0p025_seed2026101603_UNLABELED.extxyz
- Input SHA-256: da3820c477023e0bb423f87550ed82ae0c2a571f53d294c055f78d7d052d5e8b
- Approved H config SHA-256: ca6f365e8b0b8ec847e329d0312ebb51381b21e996050949494f2d3e0669dfe7
- Runner SHA-256: 0d9b97d2ed138985fc75e292485288ef7755e63399e2e4ecdaaeb30804c47e8f
- Fixed 48 atoms Ti24Si8C16, fixed cell, PBE/PW500, integer standard Monkhorst-Pack (10,10,4)
  with no explicit Gamma override. Fermi 0.10 eV, spin unpolarized,
  PPCG niter=2, Mixer(.05,8,100), energy/density/eigenstates 1e-5,
  maxiter 160 and four MPI ranks.

## H offline acceptance (NO GPAW needed)

From a fresh checkout of branch chatgpt-tsc48-segmented-runner-20261016:

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
The three included unittest cases use fake checkpoint bytes, NOT DFT.

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

RUNNER_SOURCE_CREATED=YES
PYTHON_COMPILE_PASS=YES
OFFLINE_TESTS_3_PASS=YES
CHECKPOINT_WRITE_READY=UNVERIFIED
REAL_GPAW_RESTART_VALIDATED=NO
PERSISTENT_STORAGE_READY=NO
RESOURCE_READY=NO
DFT_EXECUTION_READY=NO
DFT_STARTED=NO

The code itself does not create a persistent volume, certify wallclock
limits, perform a physical DFT restart test, or launch an unattended job.
