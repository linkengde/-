# P1G long-run output storage

Raw 64^3 snapshots and restart dumps are large. `watch_p1g_outputs.sh` compresses complete `.bin` snapshots with gzip-compatible pigz and compresses prior restart dumps only after their `.complete` marker exists. The newest checkpoint stays uncompressed for immediate recovery. `p1g_postprocess.py` reads `.bin` and `.bin.gz` snapshots transparently.

If resuming from an archived checkpoint, decompress the selected `*.dump.gz` to its original `*.dump` filename first; keep its `.state` and `.complete` sidecars beside it. Milestone dumps (`p1g_peak_pool_*.dump`, `p1g_final_solidified.dump`) remain uncompressed.

After a process has stopped, resume from the last complete checkpoint by changing to its formal case directory and running `../../resume_case_latest.sh`. It refuses to start if the same solver is already running there, restores the latest checkpoint, appends logs, and continues compression. Do not run a second copy of the case while the original is active.
