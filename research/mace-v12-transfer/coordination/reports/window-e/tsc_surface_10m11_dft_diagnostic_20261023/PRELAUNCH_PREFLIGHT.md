# Prelaunch audit: Ti3SiC2 (10-11) facet diagnostic

Date: 2026-10-23 UTC
Branch: `codex/tsc-surface-dev-prismatic-tangential-20261022`
Source audit/candidate commit: `6601f3e09a780ccfa34599e043282891e2ff2df2`

## Scope and data role

The only inputs are the 69-atom symmetric Ti-terminated `(10-11)` reference and
its paired terminal-Ti `+/-0.02 A` response from the geometry-only package.
Both descend from COD 9009647 and belong to one family. They can support a
same-source, facet-held-out diagnostic; they cannot be represented as an
independent crystallographic-source blind DEV/TEST family. No data role is
assigned, and no TEST data or evaluator-held coordinates were accessed.

## Input and environment checks

- Both input-file hashes and the runner's canonical geometry hashes match the
  manifest; each has 69 atoms, Ti33Si12C24, full PBC, and no labels.
- Python 3.12.14; isolated `/workspace/.venvs/gpaw-mpi` has GPAW 26.7.0,
  ASE 3.29.0, and gpaw-data 1.2.1. No packages were installed or upgraded.
- `GPAW_SETUP_PATH` resolves to the venv's `gpaw_data/setups` directory.
  Ti/Si/C PAW SHA-256 values match the inherited surface configuration.
- GPAW reports MPI enabled and OpenMP enabled. OpenMPI is 5.0.7; a four-rank
  `gpaw python` world-size test returned ranks 0 through 3. `OMP_NUM_THREADS=1`
  is set in the private `runtime.sh`.
- The repository's prior 100-atom surface runner is input-specific; this
  package has a separate 69-atom runner and does not modify that runner.

## Resource evidence

- CPU quota: 4 cores (`cpu.max = 400000 100000`).
- Cgroup memory limit: 34,359,738,368 bytes (32 GiB).
- Persistent `/workspace` free bytes at preflight: 16,579,280,896.
- `/tmp` free bytes: 17,170,952,192; `/tmp` is temporary and is not used as the
  result/checkpoint root.
- No active GPAW calculation was present after the initialization-only test.
- A 4-rank initialization-only attempt with the target geometry and locked
  PBE/PW500, 6x6x1 settings was stopped after 12 min 6 s before it returned an
  initialized calculator. Cgroup memory peak was 11,443,634,176 bytes. It did
  not enter SCF, produce a label, or write a checkpoint. This is a cost signal,
  not an SCF result or a completed memory benchmark.
- The largest nearby archived GPW file inspected was 1,173,989,857 bytes for a
  different 46-atom, Gamma-point diagnostic. It is not a valid direct estimate
  for this 69-atom 6x6x1 case. The exact checkpoint size therefore remains
  unknown. The runner enforces a 4 GiB post-checkpoint reserve and a 1 GiB
  transient reserve; if these cannot be maintained it stops while preserving
  the previous committed checkpoint.

## Decision

`INPUT_HASH_AND_PAWS = PASS`
`FOUR_RANK_MPI_ABI_SMOKE_TEST = PASS`
`INITIALIZATION_PREFLIGHT = INCOMPLETE_OPERATOR_STOP`
`CHECKPOINT_SIZE = NOT_MEASURED`
`INDEPENDENT_SOURCE_DEV_TEST = NOT_AVAILABLE_FROM_AUDITED_SOURCES`
`SCF_STARTED = NO` at the time this prelaunch record was written.

The fixed-geometry runner is locked to one structure per invocation. It has
memory and walltime safe-stop checks and versioned, atomic checkpoint metadata.
The reference is first; the response may start only after the reference result
has been validated and archived. Neither result changes any data role.
