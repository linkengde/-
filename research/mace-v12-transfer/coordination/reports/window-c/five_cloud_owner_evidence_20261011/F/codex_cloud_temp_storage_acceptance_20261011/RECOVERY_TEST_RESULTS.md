# Non-scientific GitHub recovery test

- Fixture: `RECOVERY_TEST.bin`, 16,415 bytes, deterministic test bytes only; no atom data or simulation output.
- Fixture source SHA-256: `24f3c1e665f4e0f74ea736d493704b7756871e238d2956761b001aff62921c6b`.
- Upload: normal fast-forward push to `origin/cloud-f-b1-dft`.
- Remote fetch: a newly initialized temporary Git repository fetched `refs/heads/cloud-f-b1-dft` from GitHub; fetched commit was `afb1ec415df43dbf6a48212723237d8a999fd00a`.
- Download: extracted `RECOVERY_TEST.bin` from the fetched commit using Git.
- Restore: copied the downloaded bytes to a second fresh path under `/tmp`.
- Verification: source, downloaded file, and restored copy were each 16,415 bytes and had the same SHA-256. `cmp` succeeded for both downloaded and restored copies.
- Result: upload, remote download, SHA-256 integrity, and byte recovery PASS for this synthetic fixture.

This does not exercise the sampler's checkpoint parser, scientific outputs, cross-configuration permissions, or repository retention policy. Those remain subject to C's signed Codex-specific release configuration.
