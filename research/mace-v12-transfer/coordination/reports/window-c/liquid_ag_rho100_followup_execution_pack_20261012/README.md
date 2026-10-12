# Ag64 rho100 follow-up execution pack

Status: **prepared for C review; execution is HOLD**. This package does not
start MD, generate atoms, evaluate the EAM calculator, or alter the archived
step1000 checkpoint. The packaged config has `execution_authorized=false`.

The package follows C's frozen strategy commit `53cc7bb8f33cf94568212b3b36a04c7f0a4a5709`
and preserves the original scientific lock and step1000 checkpoint. The source
checkpoint SHA is
`c820437a2512c4b3d7dc9e60bb60ac3ced6588aa1f2a629d58cfdabb24f34046`;
the EAM identity is checked by SHA only and remains outside Git archives.

## Frozen phase order

C's machine-readable strategy has four sequential phases. It requires an
initial 1600 K hot-liquidification stage; after its two-window gate passes,
exactly 10000 steps (5 ps) of linear cooling from 1600 K to 1250 K; then a
1250 K hold with a 2 ps discarded transient and a maximum 10 ps equilibration;
only after its two-window gate passes does it permit a 5 ps candidate window.
The user's short Stage B description combines the 1250 K transition and
candidate sampling, but the runner keeps C's explicit cooling and equilibration
phases intact. It will not skip those phases.

Every internal chunk is at most 1000 steps and 3600 wall seconds. At its last
completed step the runner writes positions, momenta, all three PCG64 states,
phase state, logs, trajectory and hashes. It then makes an immutable D-branch
archive, pushes without force, fresh-clones/fetches, restores to a new path and
compares all member hashes. Automatic continuation is allowed only after every
archive and restore check passes. The hard ceiling is global step 60000 / at
most 59000 added steps from step1000.

## Analysis implementation and remaining C gate

The module reports the requested temperature mean and standard deviation,
global Q6, CNA crystalline fraction and largest cluster, coordination, RDF,
minimum distance, stress drift, and multiple-origin MSD. Its locked numerical
gates are read from the C strategy lock. The implementation uses periodic image
offsets for neighbors; global Q6 and the FCC/HCP common-neighbor classifier use
the 3.8 Å first-shell cutoff. RDF values use 0.05 Å bins. Block bootstrap uses
2000 deterministic, circular moving-block replicates with block length derived
from the initial-positive-sequence autocorrelation estimate; MSD uses the
combined 4 ps interval for displacement and multiple origins for a long-lag
slope.

C's frozen lock also requires a “broad first peak” and “no persistent sharp
fcc/hcp multishell order,” but gives no numerical threshold for those terms.
The runner therefore requires C's execution release to include a hashed
`rdf_shape_policy` with the minimum first-peak FWHM, secondary-shell search
start, and maximum secondary-peak prominence. Until that policy is signed, the
liquid gate returns HOLD and automatic phase transitions cannot run. This
preserves the frozen physical thresholds instead of inventing a new one.

## Checks

`test_execution_pack.py` exercises numerical analysis helpers, automatic
1000-step chunk boundaries, the frozen phase transitions and terminal limits,
archive-failure continuation blocking, abnormal synthetic segment stop, a
read-only restore of the archived step1000 checkpoint, and a synthetic local
bare-Git push/fresh-fetch/restore with per-file SHA comparison. The runner uses
the same stage-transition and archive-gate functions exercised by these tests.
The tests do not construct a scientific Ag structure, attach a calculator, or
integrate a step. Package publication separately verifies the actual D GitHub
branch by fresh fetch and file SHA. A real run still requires a C release from
its stated commit,
binding the exact config, runner, criteria module, archive backend, environment,
potential, checkpoint, D branch head, and bounded plan.

Read-only preflight, after activating the pinned environment:

```sh
source /workspace/onboarding/runtime.sh
/workspace/.venvs/gpaw-mpi/bin/python followup_runner.py --preflight \
  --config FOLLOWUP_EXECUTION_CONFIG.json \
  --config-sha256 <SHA256SUMS.txt entry>
```

Do not use `--run` unless C has published the execution release and its commit,
blob path, and SHA are verified. The original checkpoint, potential, density,
cell, timestep, friction, RNG stream states and atom order are immutable inputs.
