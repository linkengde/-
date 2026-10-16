#!/usr/bin/env python3
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import age_reassemble as relay


AGE = os.environ.get("AGE_BIN", "age")
AGE_KEYGEN = os.environ.get("AGE_KEYGEN_BIN", "age-keygen")


class AgeReassembleTests(unittest.TestCase):
    def test_independent_chunk_round_trip_and_hashes(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            identity = root / "identity.txt"
            subprocess.run([AGE_KEYGEN, "-o", str(identity)], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            recipient = subprocess.check_output([AGE_KEYGEN, "-y", str(identity)], text=True).strip()
            raw = (bytes(range(256)) * 4096) + b"TSC48 synthetic relay test"
            plain = root / "plain.bin"
            encrypted = root / "part.age"
            recovered = root / "recovered.bin"
            plain.write_bytes(raw)
            subprocess.run([AGE, "-r", recipient, "-o", str(encrypted), str(plain)], check=True)
            raw_sha = hashlib.sha256(raw).hexdigest()
            cipher_sha, cipher_bytes = relay.sha256_file(encrypted)
            relay.decrypt_verified(AGE, identity, encrypted, recovered, raw_sha, len(raw))
            self.assertEqual(recovered.read_bytes(), raw)
            self.assertGreater(cipher_bytes, len(raw))
            self.assertEqual(cipher_sha, relay.sha256_file(encrypted)[0])

    def test_manifest_rejects_non_sequential_or_oversized_parts(self):
        base = {"schema": relay.SCHEMA, "repo": relay.REPO,
                "encryption": "age-v1-per-chunk", "source_bytes": 1,
                "source_sha256": "0" * 64,
                "parts": [{"index": 1, "commit": "a" * 40,
                           "path": "payload/part_000000.age",
                           "plain_bytes": 1, "plain_sha256": "0" * 64,
                           "age_bytes": 10, "age_sha256": "1" * 64}]}
        with self.assertRaisesRegex(ValueError, "PART_ORDER_INVALID"):
            relay.validate_manifest(base)
        base["parts"][0]["index"] = 0
        base["parts"][0]["plain_bytes"] = 32 * relay.MIB + 1
        base["source_bytes"] = 32 * relay.MIB + 1
        with self.assertRaisesRegex(ValueError, "PLAINTEXT_PART_EXCEEDS_32_MIB"):
            relay.validate_manifest(base)


if __name__ == "__main__":
    unittest.main(verbosity=2)
