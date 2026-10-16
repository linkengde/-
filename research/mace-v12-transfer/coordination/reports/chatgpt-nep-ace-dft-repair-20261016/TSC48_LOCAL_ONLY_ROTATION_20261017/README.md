# TSC48 H1 local-only checkpoint rotation

This branch adds an opt-in storage rotation layer around the frozen V4 GPAW
runner. It does not change the input structure or any scientific parameter.

`LOCAL_ONLY_ROTATION` accepts the user's explicit risk that loss of H1 may make
the calculation unrecoverable. The mode still refuses to remove a checkpoint
unless a newer finalized `mode=all` GPW is current, its SHA and byte count pass,
the active-run pointer has advanced, no process has the predecessor open, and
the small pre-delete progress record has been pushed to and read back from
GitHub. The predecessor is removed only after these checks. A post-delete record
captures the observed disk-space change. Any inconsistency fails closed.

At most two full GPWs exist during a save: the current recovery point and the
new `.writing.gpw`; after publication and remote metadata verification, the
predecessor is retired. GPW data itself is never uploaded. The code uses a
bounded-memory SHA-256 reader and advises the kernel to evict checksum pages.

The archive is an ordinary small-metadata branch in the existing public
repository. It contains hashes, byte counts, iteration summaries and software
provenance only, never GPW bytes, atomic coordinates or credentials. The active
checkpoint remains local to H1; this mode does not promise cross-window or
cross-session recovery.

Run the targeted offline checks with:

```sh
python3 -m unittest discover -s TSC48_LOCAL_ONLY_ROTATION_20261017/tests -v
```

The adapter must be launched only with `TSC48_LOCAL_ONLY_ROTATION=YES`,
`TSC48_ARCHIVE_REPO` set to a clean checkout of this branch, and
`TSC48_ARCHIVE_BRANCH=h-tsc48-local-only-rotation-20261017`. The normal H
resource clearance remains mandatory; local-only mode does not suppress the
disk or wall-clock gates.
