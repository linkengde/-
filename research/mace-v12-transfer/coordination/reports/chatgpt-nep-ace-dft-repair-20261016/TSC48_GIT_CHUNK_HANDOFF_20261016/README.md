# TSC48 Git chunk handoff — offline tool, NOT execution approval

The tool `chunk_transfer.py` provides **offline** splitting, SHA-256 verification,
and atomic reconstruction of a GPW or synthetic file. It does NOT push to GitHub,
run GPAW, edit the V4 runner, or remove checkpoint files. No real GPW currently
exists to transfer. The source scientific lock remains the frozen V4 commit
`5b67aa3002b1e3c197d7d511e9e075586c942095`.

## Offline 5 MiB demonstration

```bash
python -m unittest -v test_chunk_transfer.py
python -c 'from pathlib import Path; Path("sample.bin").write_bytes(bytes(range(256)) * (5 * 1024 * 1024 // 256))'
python chunk_transfer.py split sample.bin sample_parts --chunk-mib 1
python chunk_transfer.py verify sample_parts/manifest.json
python chunk_transfer.py join sample_parts/manifest.json restored.bin
sha256sum sample.bin restored.bin
```

For an actual GPW, default chunk size 32 MiB. The GitHub normal-Git hard limit
per file is 100 MiB, but uploads over 50 MiB trigger a warning. A rough
11.5 GiB GPW will produce ~368 x 32 MiB parts. GitHub limits push size to
2 GiB, so a normal Git upload must use small batches of commits and pushes,
not one giant commit/push. ***Normal Git is not a good long-term storage backend
for recurring multi-GB GPW versions***. GitHub recommends repos <1GB ideal,
<5GB strongly. Removing committed parts in later commits does not eliminate
blobs from Git history; a repeated handoff could quickly grow to tens of GB.

**Privacy:** As independently verified, `linkengde/-` is currently PUBLIC.
Do not upload unpublished GPWs to this repository without an explicit decision
to make these files public. `git push --dry-run` alone does not prove actual push.

**Required before production:** independent H-to-GitHub actual push and new-H
actual download with 5 MiB synthetic file; complete original SHA, each part SHA,
manifest hash and provenance bound to approved input/config/V4 runner;
measured workspace peak at source and destination; and independent GPAW warm
restart trial. The frozen V4 A/B staged writer can require THREE complete
GPW copies; simply splitting files or opening additional H windows does not
fix that >28.83 GB single-window disk requirement. A safe alternative would
need a separately reviewed remote-backed checkpoint/slot policy.

This program intentionally never deletes original files or removes chunks.
It also does not include any Git credentials, upload automation, or external
network access. Credentials belong in the H environment, not the repository.
