import hashlib
from pathlib import Path
import tempfile
import unittest
from chunk_transfer import split_file, load_manifest, verify_parts, join_file


class ChunkTransferTests(unittest.TestCase):
    def test_roundtrip_multi_piece(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "checkpoint.gpw"
            source.write_bytes((bytes(range(256)) * 12345)[:3 * 1024 * 1024 - 100])
            data = split_file(source, root / "chunks", chunk_mib=1)
            self.assertEqual(len(data["parts"]), 3)
            self.assertEqual(verify_parts(root / "chunks" / "manifest.json"), data)
            answer = join_file(root / "chunks" / "manifest.json", root / "out.gpw")
            self.assertEqual(answer["sha256"], hashlib.sha256(source.read_bytes()).hexdigest())
            self.assertEqual(source.read_bytes(), (root / "out.gpw").read_bytes())

    def test_tamper_and_no_partial_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "checkpoint.gpw"
            source.write_bytes(b"secret-data" * 100)
            split_file(source, root / "chunks", 1)
            p = root / "chunks" / "part_000000.bin"
            p.write_bytes(b"tampered-data")
            with self.assertRaisesRegex(ValueError, "PART_SIZE_MISMATCH"):
                verify_parts(root / "chunks" / "manifest.json")
            with self.assertRaises(ValueError):
                join_file(root / "chunks" / "manifest.json", root / "out.gpw")
            self.assertFalse((root / "out.gpw").exists())
            self.assertFalse((root / "out.gpw.writing").exists())

    def test_invalid_manifests_and_no_overwrites(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "in.gpw"
            source.write_bytes(b"abc")
            split_file(source, root / "chunks", 1)
            (root / "existing.gpw").write_bytes(b"untouched")
            with self.assertRaisesRegex(ValueError, "DESTINATION_ALREADY_EXISTS"):
                join_file(root / "chunks" / "manifest.json", root / "existing.gpw")
            self.assertEqual((root / "existing.gpw").read_bytes(), b"untouched")
            with self.assertRaisesRegex(ValueError, "OUTPUT_DIRECTORY_ALREADY_EXISTS"):
                split_file(source, root / "chunks", 1)
            data = load_manifest(root / "chunks" / "manifest.json")
            data["parts"][0]["file"] = "../../oops"
            import json
            (root / "chunks" / "manifest.json").write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "BAD_PART_NAME"):
                verify_parts(root / "chunks" / "manifest.json")


if __name__ == "__main__":
    unittest.main()
