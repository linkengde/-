# Ag64 rho100 authorization interpretation update

**Status:** `ACCEPTED`
**Purpose:** Clarify output-file naming and the already-signed `PENDING` / `QA_HOLD` continuation behavior.
**This is an interpretation supplement only. It is not a new execution authorization.**

## Bound records

- D production commit: `e228f4624e4b022023b995f12b7e05b28c92cd94`
- D runner SHA-256: `f26be9819baedc6cad5c3d1ed36ee0e7dba33be2786b09ce8bacf826e3a1c6bd`
- QA criteria SHA-256: `dadf5647cd1ffb68fd421e8f2742117c604b324e20443c0e2551e7208ae21c69`
- Production config SHA-256: `3af1d791a0e8afb601bfd10822a56ad09b2fe0d87580e5e16758057674c86c1c`
- Original C authorization commit: `e5a4d6c66af12d4eedce85644487cbed61a04eee`
- Original authorization file SHA-256: `c29842f9bfeb1c2e2013f3a077ba2bb6492a996afb00d0d3c06b77384d995b30`

## 1. Output file names

The launch-note names `temperature.csv` and `pressure.csv` are descriptive labels, not required literal filenames. The canonical files emitted by the signed D runner are:

- `thermal.csv`: temperature, target temperature, and energy records. The runner records these with global step and time.
- `stress.csv`: all six stress components, likewise keyed by global step and time.

The QA code joins the two canonical files by global step. It derives pressure from the diagonal stress components using the existing locked convention, `P = -(sigma_xx + sigma_yy + sigma_zz) / 3`. A separate `pressure.csv` is not required. Do not rename the files or add duplicate CSVs for this clarification.

**Decision:** the mapping `temperature.csv` → `thermal.csv` and `pressure.csv` → `stress.csv` plus derived pressure is accepted. No runner or scientific-file change is needed.

## 2. QA status and continuation

- `PENDING`: the locked QA does not yet have enough eligible data or statistical support to decide. It is not a pass and does not trigger a stage transition. The signed runner leaves the current stage unchanged and may proceed to the next authorized chunk only after that chunk has been checkpointed, archived, fresh-fetched, hash-verified, and restored successfully, and while all existing step, wall-time, resource, and stage limits remain satisfied.
- `QA_HOLD` (the runner's explicit HOLD gate): required input or a metric is missing, malformed, or otherwise fails its locked QA validity checks. The runner changes to `HOLD_QA_REQUIRED_INPUT_OR_METRIC` and stops; it must not begin another chunk.
- Other statuses, including `PASS` and `FAIL_WINDOW`, retain their existing signed runner and criteria behavior. This supplement does not redefine them.

**Decision:** `PENDING` continues within the already-authorized bounded run; `QA_HOLD` stops. `PENDING` must never be reported as `PASS` or as a confirmed liquid candidate.

## Scope

This supplement only resolves the two interface interpretations above. It does not modify or replace the original authorization, D commit, runner, criteria, checkpoint, EAM potential, configuration, output data, scientific parameters, QA thresholds, or authorized execution bounds. The original signed authorization remains controlling; all its existing gates and limits remain in force.

Code basis reviewed at D commit `e228f4624e4b022023b995f12b7e05b28c92cd94`: the production runner reads `thermal.csv` and `stress.csv`, returns an ordinary chunk-complete status for `PENDING`, and converts the explicit `QA_HOLD` gate into a terminal hold state before continuation is considered.
