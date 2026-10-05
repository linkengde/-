# Provenance-aware V15 builder revision — Window B

This is an isolated report artifact. A decides production integration and whether to authorize any numerical screening. Nothing here authorizes training. Dependencies: existing ASE and NumPy; no installation, DFT or model imports.

`build_v15_draft.py` is a copy of the published B builder, retaining parent split hashes, four LCAO exclusions, three residual archive/hash/convergence checks, finite labels, exact serialization checks, historical valid/test byte copying, composition rank, geometry leakage diagnostics, and frozen blind INPUT checks. New features: separately stored `source_method`/`dft_method`, hash-bound trace evidence, evidence tiers, and explicit sole-unresolved-row exclusion. Missing settings remain unknown. REF energy is never replaced by free energy.

Run from `/tmp` (commands actually exercised):

```bash
cd /tmp
/workspace/.venvs/gpaw-mpi/bin/python -B /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v15_provenance_aware_builder_revision/build_v15_draft.py \
  --repo-root /workspace/- \
  --output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v15_provenance_aware_builder_revision/baseline_30
/workspace/.venvs/gpaw-mpi/bin/python -B /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v15_provenance_aware_builder_revision/build_v15_draft.py \
  --repo-root /workspace/- \
  --output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v15_provenance_aware_builder_revision/exclude_unresolved_29 \
  --exclude-unresolved-agti
```

Both output directories are now nonempty; rerunning these exact commands intentionally refuses overwrite. For a new run choose a new empty child of this report directory. Production paths and nonempty destinations are rejected. Builder failures retain `diagnostics_failure.json` in the selected isolated output directory; do not suppress failures or overwrite diagnostics.

The exclusion selector requires the config `AgTi_2p3871A_periodic_PW_PBE_force_only`, V14 source frame 13, source SHA256 `3323cdbd3081fc8d91673662c5f83c8da93ba3b4d7a23428f6b70d0e7a4c9faf`, and trace source path to match exactly once. It does not exclude the other 14 partially documented rows. Four historical force-only LCAO rows remain excluded in both scenarios.

Independent scenario validation (reads only the prior B draft, these scenario files and explicit hashed sources; no blind label files):

```bash
cd /workspace/-
/workspace/.venvs/gpaw-mpi/bin/python -B research/mace-v12-transfer/coordination/reports/window-b/v15_provenance_aware_builder_revision/verify_scenarios.py
sha256sum -c research/mace-v12-transfer/coordination/reports/window-b/v15_provenance_aware_builder_revision/SHA256SUMS.txt
```

Each scenario contains train/valid/test structures with inherited REF labels, its own manifest, verification report and checksum list. The small split files are needed for scenario verification; no models, checkpoints or original historical datasets were recopied beyond these required outputs. `comparison.json` and `comparison_report.md` summarize the result and limitations. Root checksums cover published deliverables and exact current source inputs; check from repository root. Scenario checksums use filenames relative to their scenario directory.
