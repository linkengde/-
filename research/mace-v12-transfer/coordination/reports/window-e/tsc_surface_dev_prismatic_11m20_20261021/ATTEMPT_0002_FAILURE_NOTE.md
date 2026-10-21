# Attempt 0002 failure note

Attempt `0002` stopped at the scheduled checkpoint callback for logical SCF iteration 20. The SCF log contains iterations 1–19 and no convergence marker. The checkpoint callback raised `KeyError: 'disk_stop_guard_bytes'` because the resource guard was present in supervisor preflight but missing from `SCIENTIFIC_CONFIG.json`.

No checkpoint or DFT labels were written. The attempt directory is retained unchanged. This was a runner configuration defect, not evidence of insufficient disk or memory; the supervisor measured about 16.95 GB peak cgroup memory and about 16.64 GB free workspace at exit.

The guard was added as 4,294,967,296 bytes, matching the pre-existing supervisor hard disk guard. Startup validation was added so a missing or inconsistent value fails before SCF. The input structure, coordinates, cell, PBC, PAW files, k-point mesh, functional, cutoff, occupations, eigensolver, mixer, and convergence thresholds remain unchanged. Python compilation and fresh resource preflight passed after the patch.

Attempt `0003` is a new calculation from the original frozen input because attempt `0002` produced no verified checkpoint to resume. No files from attempt `0002` were overwritten.
