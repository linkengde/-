# V17 holdout DFT execution pack (preparation only)

This package prepares nine geometry-only PW-PBE single-point inputs: three proposed development-validation points and six proposed withheld-test points. **It has not launched DFT.** Every label remains `blocked_pending_A_role_review_and_owner_assignment`; no authorization file is included or self-issued.

## Proposed roles and owners

| Role | Proposed owner | Labels |
|---|---|---|
| Development validation | A | `AgC_v17_dev_validation_01`, `AgSi_v17_dev_validation_01`, `AgTi_v17_dev_validation_01` |
| Withheld test | B | `AgC_v17_test_distance_01`, `AgC_v17_test_transverse_framework_01`, `AgSi_v17_test_distance_01`, `AgSi_v17_test_transverse_framework_01`, `AgTi_v17_test_distance_01`, `AgTi_v17_test_transverse_framework_01` |

This split is only a proposal. A must review and accept every role, assign each owner, and issue a separate authorization JSON before any future launch. The preflight refuses to proceed if that file is absent, has the wrong role-manifest/input hashes, or lacks A's role and owner decisions. Do not create an authorization file in this B package.

## Geometry and screening provenance

The nine `inputs/*.extxyz` files are byte-for-byte copies of the frozen candidate geometries in `v17_split_geometry_design/final_set/geometries/`. `input_manifest.json` records their SHA256 hashes, ordered atom IDs/elements, formula, cell/PBC, canonical geometry identity, parent-input hashes, role, proposed owner, marked pair, species-specific minimum distances, and per-candidate duplicate/contact screening. `source_inventory.json` records geometry-only source hashes. The historical data-set hashes listed as lineage are provenance references from the earlier geometry audit, not inputs to this pack.

Preparation deviation: an early broad recursive search used to locate the pinned source settings also matched some existing V17 calculation log/summary text. I stopped that search and did not copy or use those result values; no further result folders were read. This pack does not analyze those results. See `readiness_report.md` for the disclosure.

The prior geometry audit reported zero exact/near duplicates against its inventory and pairwise-separated proposed roles. The candidates still inherit existing small-cluster motifs. They do not establish new morphology or thermal-state coverage.

## Pinned calculation recipe

If A later authorizes execution, the drafts preserve fixed-geometry single-point GPAW 26.7.0, PBE, plane-wave cutoff 500 eV, Gamma-only `(1,1,1)`, Fermi-Dirac smearing 0.1 eV, Pulay mixer (`beta=0.05`, `nmaxold=8`, `weight=100`), energy/density/eigenstate thresholds `1e-5`, four MPI ranks, and one OpenMP/BLAS thread per rank. The native extrapolated energy is the primary energy and free energy is stored separately when GPAW provides it. `state.gpw` remains in the local run root and is never copied to the compact archive.

## Future use after explicit A authorization

Only after A publishes the accepted role manifest and a separate authorization file may the assigned owner run the guarded draft. The authorization must include `status: authorized`, `reviewed_by: window-a`, the exact role-manifest hash, exact owner/role/input-hash maps, and (for every withheld-test label) the frozen selected-model hash, selection-record hash and freeze timestamp before A reads test labels. Any role or hash mismatch stops before MPI/GPAW is started.

The eventual runner requires explicit `--workspace-root`, `--pack-root`, `--input-root`, `--run-root`, `--archive-root`, `--authorization`, `--owner`, and `--labels` arguments. It refuses existing per-label run/archive directories, stale staging paths, and overlapping MPI/DFT processes; it checks four available cores, GPAW MPI/ScaLAPACK support and one thread per rank. It does not claim tasks or modify A's queue.

Example shape only; the authorization path is intentionally not present in this package:

```bash
bash run_holdout_labels.sh \
  --workspace-root /workspace/- \
  --pack-root /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v17_holdout_dft_execution_pack \
  --input-root /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v17_holdout_dft_execution_pack/inputs \
  --run-root /tmp/v17_holdout_run \
  --archive-root /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v17_holdout_dft_execution_pack/archives \
  --authorization /path/to/a-issued-authorization.json \
  --owner window-b \
  --labels AgC_v17_test_distance_01 AgC_v17_test_transverse_framework_01
```

For withheld tests, B's allowed responsibility is only an explicitly assigned DFT run and mechanical archive-integrity verification. B must not analyze test energies/forces or run model inference/scoring. A records the selected checkpoint before reading the withheld labels and owns interpretation. These labels must not enter model or parameter selection.

## Package files

- `input_manifest.json`: roles, proposed owners, methods, input hashes and geometry identity.
- `source_inventory.json`: candidate and geometry-parent hashes; no DFT outputs.
- `prepare_pack.py`: geometry-only manifest builder used for this preparation.
- `execution_preflight.py`: future role/hash/geometry/authorization guard; does not import GPAW.
- `run_holdout_labels.sh`: guarded future four-rank runner draft.
- `dft_single_point.py`: future GPAW single-point draft; GPAW is imported only after explicit authorization and geometry checks.
- `archive_verify.py`: future convergence/finite-value/identity verifier and compact archive writer. For withheld tests it reports pass/fail only and prints no label values.
- `readiness_report.md`: completed preparation checks and limits.
- `SHA256SUMS.txt`: package hashes, excluding itself.
