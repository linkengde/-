import io
from pathlib import Path
import tempfile
import unittest
from chunk_transfer import split_file
from stream_join_github import stream_join


class RemoteTests(unittest.TestCase):
    def test_streaming_no_chunk_copies(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp)
            original=folder/'in.gpw'
            original.write_bytes(bytes(range(256)) * 9000)
            split_file(original, folder/'chunks', 1)
            def opener(url):
                return io.BytesIO((folder/'chunks'/url.rsplit('/', 1)[-1]).read_bytes())
            dest=folder/'out.gpw'
            result=stream_join(folder/'chunks'/'manifest.json',
                'https://raw.githubusercontent.com/linkengde/-/abc123/test/', dest, opener=opener)
            self.assertEqual(result['status'], 'SHA_MATCH')
            self.assertEqual(dest.read_bytes(), original.read_bytes())
            self.assertFalse((folder/'out.gpw.writing').exists())

    def test_remote_tamper_fail_and_no_output(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp)
            (folder/'in.gpw').write_bytes(b'x'*300)
            split_file(folder/'in.gpw',folder/'chunks',1)
            dest=folder/'out.gpw'
            with self.assertRaisesRegex(ValueError, 'REMOTE_PART_HASH_MISMATCH'):
                stream_join(folder/'chunks'/'manifest.json',
                    'https://raw.githubusercontent.com/linkengde/-/abc123/test/', dest,
                    opener=lambda u: io.BytesIO(b'z'*300))
            self.assertFalse(dest.exists())
            self.assertFalse((folder/'out.gpw.writing').exists())

    def test_https_only_and_existing_destination(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp)
            (folder/'in.gpw').write_bytes(b'abc')
            split_file(folder/'in.gpw',folder/'chunks',1)
            with self.assertRaisesRegex(ValueError,'ONLY_HTTPS_RAW_GITHUB_SUPPORTED'):
                stream_join(folder/'chunks'/'manifest.json','http://example.com/evil/',folder/'new')
            (folder/'out').write_bytes(b'no-clobber')
            with self.assertRaisesRegex(ValueError,'DESTINATION_ALREADY_EXISTS'):
                stream_join(folder/'chunks'/'manifest.json',
                    'https://raw.githubusercontent.com/linkengde/-/abc123/test/',folder/'out')
            self.assertEqual((folder/'out').read_bytes(),b'no-clobber')

if __name__=='__main__':
    unittest.main()
