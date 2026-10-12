# Ag64 rho100 bounded continuation runner, QA v2.1

This D-owned package binds the automatic continuation runner to C's QA v2.1 science lock (`848f8db72d83cea91dcf4cd42ee9fbbf49a80e7b5e3759cb4f535be8cc5873d2`). It keeps the existing step-1000 checkpoint and WMH2006 EAM file read-only.

The historical 11-frame regression compares the legacy coherent/global Q6 implementation to `0.08871920638262233 ± 0.01`. The liquid screen requires global Q6 `< 0.20` in at least 80% of frames in each of the two final adjacent 2 ps windows. Mean local Q6 is reported only; it has no threshold or trend gate. PTM FCC/HCP at RMSD 0.10, periodic cluster limits, finite-cell RDF and coordination checks, MSD diffusion, and thermal stationarity remain required. Insufficient or unresolved required data returns PENDING/HOLD.

The runner has a preparation-only configuration with `execution_authorized=false`. It performs live, read-only preflight and will refuse to run unless C publishes a separate authorization blob in a specified C commit. That blob must match the exact D branch head, config, runner, QA lock, environment, checkpoint, potential, and step/time caps. The runner does not treat the QA science lock itself as execution authorization.

The authorized path is bounded to the existing trajectory: step 1000 onward, at most step 60000, no more than 1000 added steps or 3600 wall seconds per chunk. Before each chunk it checks the D branch against GitHub. At each boundary it writes a checkpoint, uploads an immutable segment, fresh-fetches and verifies every SHA, reloads the checkpoint positions/momenta/RNG/stage state, and only then allows the next chunk. No new source structure, RNG, velocity, or potential is generated.

Run the isolated tests with:

```sh
/workspace/.venvs/gpaw-mpi/bin/python -m unittest -v test_production_runner_v2_1.py
```

Run the read-only live preflight with the SHA recorded in `FOLLOWUP_EXECUTION_CONFIG_V2_1.json`:

```sh
/workspace/.venvs/gpaw-mpi/bin/python followup_runner_v2_1.py \
  --config FOLLOWUP_EXECUTION_CONFIG_V2_1.json \
  --config-sha256 <config SHA256> \
  --preflight
```

No formal MD has been run as part of this package. A separate real-GitHub integration harness uploads only fabricated test files; it does not include an EAM potential or scientific checkpoint.
