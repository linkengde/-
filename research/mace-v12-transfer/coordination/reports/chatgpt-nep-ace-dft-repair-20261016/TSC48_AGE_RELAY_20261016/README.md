# TSC48 age-encrypted checkpoint relay

This directory contains a receiver that downloads one independently encrypted
age chunk at a time from immutable commits in the public `linkengde/-` repository,
verifies the ciphertext SHA-256, decrypts with an identity local to H2, verifies
the plaintext part SHA-256, and atomically publishes the reconstructed GPW only
after the full-file byte count and SHA-256 match the manifest.

The receiver does not read, print, upload, or copy the age identity. Each plaintext
part is at most 32 MiB; the temporary encrypted and plaintext parts are removed
after each part has been checked. The output `.writing` file is not a valid GPW
and is never renamed to the destination unless every check passes.

## Synthetic H2 key-pair check

`test_manifest.json` refers to one 1 MiB synthetic payload encrypted to the H2
recipient public key. H2 should fetch this manifest and run the receiver before
the real GPW upload begins:

```sh
curl -fL '<manifest raw URL from the handoff message>' -o test_manifest.json
python3 age_reassemble.py test_manifest.json \
  --identity '<H2-local age identity path>' \
  --output /tmp/age-key-pair-check.bin
```

Expected result SHA-256:
`fbbab289f7f94b25736c58be46a994c441fd02552cc6022352e3d86d2fab7c83`

H2 must report the actual output SHA and confirm that its private key decrypted
the test. H1 has only the public key and cannot attest to private-key matching.

## Manifest format

The JSON schema is `h2-age-gpw-relay-v1`. `source_bytes` and `source_sha256`
refer to the original plaintext GPW. Every ordered `parts[]` record contains
`index`, immutable `commit`, repository-relative `path`, `plain_bytes`,
`plain_sha256`, `age_bytes`, and `age_sha256`. Ciphertext locations use raw HTTPS
URLs derived from the fixed repository, commit, and path.

## Source provenance

The branch is based on the previously accepted relay commit
`43693387544e9a1e50169b1ae92918438b3ef144`. V4 remains frozen. The transferred
checkpoint's `current.json` records the executed adapter SHA
`fc32f931dc9ac336ee1ab11057293ae7df013b8be781a2de4d4d30a8e4a15dbb`; that
adapter was an uncommitted local revision and differs from the adapter in the
accepted commit (`53553bcc91e6e8cde86788a0c4588cd92e07cfe7d4ae705407d1e6bb7b46ba77`).
`executed_relay_runner.py` preserves the exact adapter bytes that generated the
checkpoint. Its inclusion does not change the scientific parameters.

## Local tests

The test suite generates a disposable local age key in a temporary directory,
round-trips a synthetic chunk, checks the hashes, and rejects invalid manifests.
It does not use the H2 identity, run GPAW, or create scientific data.
