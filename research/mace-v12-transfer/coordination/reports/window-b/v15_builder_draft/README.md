# Isolated V15 builder draft (Window B)

A reviews and integrates this draft. The generated files are report artifacts, not a production dataset or authorization to train. No DFT, MACE inference/training, MD or TTM was run. V15 blind label files were never read.

`build_v15_draft.py` requires explicit repository-root and output-directory arguments, independent of its own location and the caller's working directory. It permits output only in a child directory of this B report directory and refuses a nonempty output directory. It checks the V14 split hashes and counts, excludes the four documented force-only LCAO records, validates the three approved residual archives, and copies historical valid/test bytes verbatim.

It records actual energy/force fields rather than inferring them from historical config tags. Unknown inherited methods remain unknown. No wide-span unlabeled proposal is included. See `verification_report.md`, `exercise_checks.json`, and the generated manifest for limitations and checks.

This exact command succeeded from `/tmp`:

```bash
cd /tmp
/workspace/.venvs/gpaw-mpi/bin/python \
  /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v15_builder_draft/build_v15_draft.py \
  --repo-root /workspace/- \
  --output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v15_builder_draft/generated_data
```

The existing `generated_data` directory is nonempty and the command now intentionally refuses to overwrite it. For a review rerun, choose a new empty child such as `review_run_02` with the same two explicit arguments. Dependencies are the already installed ASE and numpy; no DFT/MACE imports or installations are needed. A production path is intentionally rejected.

Outputs:

- `generated_data/train.extxyz`, `valid.extxyz`, `test.extxyz`: 30/2/3 draft frames.
- `generated_data/dataset_manifest.json`: 35 per-frame provenance records, four exclusions, source/output hashes, counts and composition matrix.
- `generated_data/verification_report.json`: archive checks, blind-input isolation, historical and retained split overlaps.
- `exercise_checks.json`: independent parsing/hash/count checks and refusal/geometry guard checks.
- `SHA256SUMS.txt`: all deliverables plus every source read, excluding itself.

From `/workspace/-`, verify the publication with:

```bash
sha256sum -c research/mace-v12-transfer/coordination/reports/window-b/v15_builder_draft/SHA256SUMS.txt
```

The generated-data checksum file separately uses paths relative to `generated_data`; verify it from that directory. If an archive/source validation fails during a new run, the builder exits nonzero and retains `diagnostics_failure.json` and any partial output in the isolated output directory. Investigate the blocker; do not overwrite the nonempty run directory or suppress checks.
