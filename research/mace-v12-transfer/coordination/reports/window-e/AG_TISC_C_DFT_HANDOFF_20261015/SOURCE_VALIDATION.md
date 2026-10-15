# Source validation record

All paths below were checked in the local repository before copying. The parent audit commit remained `a0e94666eb6a7a355b9c1e5bcaba58c9868cd8c2` on `cloud-e-b1-dft`; this packaging operation did not commit or change that revision.

- Audit folder: its five payload reports passed its original `SHA256SUMS.txt` check.
- 135-atom model archive: the source model manifest passed for all 22 listed entries; input and optimized structure identities agree with `RESULT.json` and the audit's stated hashes.
- Component force attribution: its original source `SHA256SUMS.txt` verified all 53 listed files from the correct `minimal_135atom_model` base directory. The two large GPW checkpoints were validated at source but intentionally not copied; component output JSON and forces retain the checkpoint SHA identifiers.
- V19 reference bundle: embedded SHA verification and `python3 verify_bundle.py` passed: 132 files and 16 verified reference cases.
- candidate02 DFT diagnosis archive: all 5 listed files passed its source SHA manifest.
- candidate02 decomposition archive: all 7 listed files passed its source SHA manifest.

Each copied file's source path, source-side SHA256 and byte count is recorded in `SOURCE_INVENTORY.csv`. The handoff-root `SHA256SUMS.txt` independently binds the exact packaged bytes. Failed/incomplete 2×2×1 data are classified as failure evidence, not labels. No missing structure is substituted or generated in this package.
