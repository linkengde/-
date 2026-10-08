# Overnight A/B controlled width workflow

Use the current checkout and registered new instance. Do not create a worktree or claim another owner's job. Existing 0.05 eV jobs are allowed to finish. The driver recognizes active queues/calculators and verified archives; it does not repeat a finished job or overwrite an unfinished run.

For new B, from `/workspace/-`:

```bash
git pull --ff-only origin main
mkdir -p /workspace/.setup
nohup python3 -u research/mace-v12-transfer/coordination/run_overnight_width_scan.py window-b --hours 9 > /workspace/.setup/v19-overnight-window-b.log 2>&1 < /dev/null &
```

Run the entry once; a local driver lock rejects duplicates. No terminal must remain open. A uses `window-a`. Each owner advances from its own verified 0.05 result to its own registered 0.20 job without waiting for the other owner. A's previous completion controller was replaced while its live DFT queue and MPI ranks were preserved.

1. Verify ownership and synchronize ordinary `main` without reset/force push.
2. Wait for the existing 0.05 job; start it only if neither run nor archive exists. Verify the complete archive and mark its task label complete.
3. Run the registered same-input/method/mesh 0.20 eV job, one per owner/machine. Queue and calculator require the verified 0.05 predecessor. Preserve four MPI ranks and one thread/rank, input hashes, method and disk guards.
4. Verify/publish compact results; retain local `state.gpw` without uploading it.
5. Compare 0.20 against both 0.10 and 0.05 at the owner's fixed mesh, then compare 5/6 at 0.20 once the other owner publishes. Publish machine-readable comparisons and report.

Inspect `/workspace/.setup/v19-overnight-window-b.json` and `.log` for B's current stage. Failure stops with the same run/checkpoint preserved; diagnose and resume, do not duplicate jobs or change gates. The nine-hour bound controls waiting and new launches; it never kills an active calculation at the deadline. A long job can finish beyond the waiting window.

These are numerical sensitivity jobs, not model validation labels. Wider smearing does not prove improved physical accuracy. No old grid rerun, sealed-label reads, dataset integration, model fitting, LAMMPS/QE, or automatic vacuum job. Review the combined width/mesh evidence before preparing the next physical-setting check.
