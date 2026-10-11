# Checkpoint identity and hash scope

The archived field called `scientific_parameters_sha256` is retained as historical provenance. It is not a digest of only physical inputs and must not be reused as the new executor identity.

The old runner computes it as SHA-256 of the canonical JSON encoding of `source_manifest.json`'s `fixed_identity` object: `json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()`. The digest was independently recomputed as `8f792d6206c1b7bae7b900e7ea81445d6cd3e837de2b5979c88bb977650ab7b2`, matching both the archived manifest and checkpoint state.

That object includes the potential metadata and file SHA; sampling settings (excluding only `max_wall_s`, `warmup_max_wall_s`, `max_rss_bytes`, `min_free_bytes`, `trajectory_archive_max_bytes`, and `checkpoint_reserve_bytes`); task row; source root; task-matrix SHA; package versions; old executor and parser source hashes; sampler root; ASE Langevin version marker; hostname; Python version; thread environment; and thermostat description. Thus this historical digest covers code and host details as well as scientific settings. In particular, its `run_sampling_codex.py` hash is `6b54d094c26bb9d4a5e4eb5997a31955eb8829a2975042e23fe10527cb943583`, independently matching the archived source file.

The new adapter identity is a separate SHA-256 of `ase_segment_adapter.py`. The test code has its own SHA-256. The unsigned candidate lock binds the new code hash, old source archive/state/manifest hashes, old runner hash, fixed step boundary, potential SHA, and exact segment settings. A future C release must publish a fresh lock SHA and adapter SHA together; changing either code or lock invalidates the prior candidate.

No hash is presented as a C signature in this package. The current lock has `execution_authorized: false`, `remote_archive_verified: false`, and `release_status: PREPARED_UNSIGNED_NOT_EXECUTABLE`.
