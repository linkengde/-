# Ag64 step2000→3000 unsigned candidate config review

## Outcome

Created `UNSIGNED_STEP2000_3000_CONFIG.json` as a new file. Existing configuration files were not modified. This is an unsigned review candidate and cannot authorize or start execution.

**CONFIG_SHA256:** `1aa03f05ed2a3662985269c2b0148542d79202f771aff690d67894e25f7c1f96`

The hash is SHA-256 over the exact UTF-8 bytes of the JSON file, including its final newline.

## Bound values

| Item | Value |
| --- | --- |
| Task | `C-D-AGL-1250K-rho097-seed01` |
| Segment | `start_step=2000`, `end_step=3000`, `max_added_steps=1000` |
| Step2000 checkpoint state SHA-256 | `5d03864f07a98f177c3c0cbdcb850158ed89a64fa2a6ab7f879e8c7ab003be2a` |
| Source archive SHA-256 | `c3a88a759e0bd258e74632aba0a68e4f1a765fcbe8324034d5bf127b5aa69965` |
| Adapter SHA-256 | `4dfc0485f90d67cb56ab9ee8c8721ee2bd749aa13027d15d6c0f99c7d2a92fbb` |
| Core SHA-256 | `adcc939fcb83b3ab3e4ed6e121e4ce3615bef73e1551c78a09f347d95985e3af` |
| Archive backend SHA-256 | `0a8d8b71ee938e768f47bc9e819c5953f7f7bb0aa6415f248916b4a300a3928f` |
| Potential SHA-256 pin | `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447` |
| Environment SHA | `14084cfa4469b19d34c7dc51b03a87971c7fb69956eab539e7a127fe69e4f6d4` |
| Environment validation branch commit | `fff2d118d2fd62d52809bfd3da991edaf51125e2` |

## Review checks

- JSON parses successfully.
- Boundary fields are exactly 2000→3000 with 1000 added steps.
- The candidate records the reported ASE 3.29.0 dry-run evidence and supplied `ENVIRONMENT_SHA`.
- Adapter, core, and archive backend SHA-256 values match the source files on `ag64-step2000-extension-dev`.
- The potential appears only as a SHA-256 pin; no potential file was loaded.
- `execution_authorized` is `false`; the release status is unsigned; no C signature is included.
- No MD, force calculation, energy calculation, or dynamics step was run for this task.

The candidate intentionally does not satisfy the signed-release gate in the adapter. C review/signature must produce a separately signed release configuration before any execution can be authorized.
