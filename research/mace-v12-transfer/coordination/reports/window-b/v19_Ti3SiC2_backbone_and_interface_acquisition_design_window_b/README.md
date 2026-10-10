# Ti₃SiC₂ backbone and Ag-interface geometry proposals

## Scope

This is a geometry and source review for `v19_Ti3SiC2_backbone_and_interface_acquisition_design_window_b`. It produces four unlabelled candidate geometries:

- Four fixed-cell Ti₃SiC₂ points: ±0.02 Å along z for one representative Ti 4f site (local ID 3) and one C 4f site (local ID 9).

No new Ag–Ti, Ag–Si or Ag–C geometry was generated. Those three interface entries are explicitly blocked in `proposal_manifest.json` and `source_novelty_ledger.csv` because the public COD CIF could not be externally verified.

Every proposal has `launch_enabled=false`, `owner=null`, and no energy or force labels. The interface method remains unset pending A's mesh, smearing, vacuum, PBC and dipole review. The Ti₃SiC₂ method hash is a reference to the existing registered 10×10×4 baseline method declaration; it does not register or authorize these four points.

No DFT/MPI, MACE inference/training, dataset integration or input-manifest changes were made. The only training-related source read was the completed V18 audit's `training_provenance.csv`; the raw training file and all development, test, holdout and sealed data were not opened.

## Source and novelty review

The Ti₃SiC₂ baseline input SHA-256 matches its registered manifest entry. The repository's seven COD parent package checksums pass, including source CIF SHA-256 `6abbf0a40208c5985cacc3b9ac3564a944702203ad57e961010a674f30d598c3`; the three existing parent geometry hashes also match the package manifest. The repository records public COD entry 9009647 and DOI `10.1016/S0022-3697(98)00226-1`.

I attempted to fetch the source from the declared COD URL twice, including an approved retry. Both requests returned HTTP 403, so the live remote CIF bytes could not be compared with the repository copy. The repository hash chain passes, but the interface source is not confirmed to the level required for new proposals. Per the task boundary, Ag–Ti, Ag–Si and Ag–C remain blocked and no derived parent or offset geometry is included. This report does **not** claim a remote byte-for-byte COD match.

Against the selected V18 train-only provenance report, the COD accession and the exact source/parent geometry hashes are absent. The interface chemistries Ag–C, Ag–Si and Ag–Ti already occur in V18 training, while the completed coverage audit found no pure Ti₃SiC₂ training frame. Some V18 rows lack explicit source hashes, so the comparison cannot prove complete structural independence. The three existing interface parents point to the **same possible COD CIF family**; they would not be three independent validation sources if later verified.

The detailed frame-derived comparison is in `source_novelty_ledger.csv`. It uses only the V18 audit report, not training labels or external archives.

## Geometry checks

The generator preserves atom order, atom IDs, symbols, cell and PBC. Each Ti/C point changes only its selected atom's z coordinate. It checks periodic species-pair distances with a conservative floor of 0.80 times the sum of ASE covalent radii. This is an overlap screen only; it does not establish stability or suitable bonding. All four bulk proposals pass the screen and are unique under an identity-preserving minimum-image check at 1×10⁻⁸ Å.

## Files and reproduction

- `inputs/*.extxyz`: four Ti₃SiC₂ geometry-only proposals; no interface geometry is included.
- `proposal_manifest.json`: source lineage, bulk operations, atom identities, method reference hash, periodic distance checks, interface blockers and disabled ownership state.
- `source_novelty_ledger.csv`: comparison against the V18 train-only provenance report.
- `validation.json`: summary of integrity and geometry checks, including the unavailable remote COD hash check.
- `prepare_proposals.py`: reproduces proposal geometries and machine-readable ledgers using ASE.
- `SHA256SUMS.txt`: hashes every deliverable in this report directory except itself.

From the repository root:

```bash
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/coordination/reports/window-b/v19_Ti3SiC2_backbone_and_interface_acquisition_design_window_b/prepare_proposals.py
sha256sum -c research/mace-v12-transfer/coordination/reports/window-b/v19_Ti3SiC2_backbone_and_interface_acquisition_design_window_b/SHA256SUMS.txt
```

The generator reads the public baseline and source manifests, hashes the three existing interface parents without transforming them, and reads the completed V18 train-only provenance report. It creates no calculator and launches no calculation.
