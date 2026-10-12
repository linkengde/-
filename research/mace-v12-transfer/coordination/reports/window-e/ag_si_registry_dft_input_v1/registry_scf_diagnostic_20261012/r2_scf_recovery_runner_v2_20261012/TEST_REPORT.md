# R2 recovery runner v2 test report

**Status:** PASS for preparation and synthetic validation. **DFT execution:** NOT AUTHORIZED and NOT STARTED.

## Executed checks

- Nine serial unit tests passed. They verify the exact Phase2 lock path/SHA, R2 source hashes, explicit `PPCG(niter=4)` parameter delivery into GPAW's new ASE calculator, the 200-step and 14,400-second limits, 40-step checkpoint policy, safe-stop checkpointing, synthetic checkpoint write failure, rejection of partial outputs, acceptance of a complete synthetic converged record, and rejection of `--run` without C authorization before output-directory creation.
- Single-rank real GPAW synthetic GPW readback passed. A synthetic H initialization (SCF loop count remained zero) wrote a GPW with `mode="all"`; GPAW `read_gpw` restored wavefunctions; SHA-256 and sidecar iteration metadata matched.
- Four-rank MPI synthetic rank-zero write-error test passed. All ranks exited the checkpoint path with a failure; no successful checkpoint or label was fabricated.
- `collector.py` was tested with partial numeric fields and no convergence marker. It returned `FAILED_SCF` and omitted total energy, free energy, forces and stress.
- Exact R2 input, structure and Phase2 lock SHA values were re-read. No tracked original Phase2 file changed.
- `launch_r2.sh` uses the GPAW MPI Python launcher, OMP=1, four MPI ranks, a 12,600-second orderly-stop threshold and the 14,400-second outer wall-clock deadline.

The tests used synthetic H/fake records only. They did not run an R2 SCF iteration, evaluate R2 energy/forces/stress, alter R1/R2 inputs, or start any scientific workload.

An optional four-rank synthetic H readback test completed once. Its final repeat remained in the synthetic LCAO initialization phase for more than one minute and had not reached the checkpoint hook; it was stopped and is not counted as a pass. The repeatable serial real-GPAW readback and four-rank write-error synchronization tests passed. This is recorded as a test limitation rather than a GPAW checkpoint failure.

## Checkpoint behavior verified

GPAW 26.7.0's `ASECalculator.hooks["scf_step"]` fires after each yielded SCF iteration. `calc.write(path, mode="all")` calls GPAW's collective GPW writer with `include_wfs=True`. The synthetic MPI round-trip verified this API and its readback. A checkpoint metadata sidecar binds its SHA-256 to the callback's completed iteration. Checkpoints are immutable per iteration; a later write does not overwrite an earlier verified file.

The runner catches a write failure at the collective boundary and broadcasts a stop to all ranks. It leaves any incomplete `.partial` file and prior valid checkpoints for review, and it never retries. The orderly wall-clock path checkpoints at the first completed iteration after 12,600 seconds and exits; GNU `timeout` is the absolute 14,400-second guard if a single iteration becomes stuck. If that outer guard fires, the latest verified periodic checkpoint remains authoritative and the collector rejects incomplete labels.

## Resource check

The current task workspace is a Codex Cloud temporary workspace. The `/workspace` filesystem reported 27,939,389,440 free bytes and 2,038,336 free inodes. Cgroup CPU quota is 4.0 cores; memory limit is 34,359,738,368 bytes, current cgroup use was 2,101,415,936 bytes, and cgroup headroom was 32,258,322,432 bytes. `/proc/meminfo` reported 33,993,160 KiB available at preflight time.

Archived R2 runtime diagnostics reported 1,024,913,408 bytes for wavefunctions and 384,342,528 bytes for PAW projectors. The checkpoint working-set estimate uses twice that combined footprint: 2,818,511,872 bytes; cgroup memory headroom exceeds it. The budget estimates each full GPW at the same 2,818,511,872 bytes. Five 40-step checkpoints plus one simultaneous temporary file, 1 GiB of logs and a 5 GiB reserve require 23,353,522,176 free bytes. The measured free-space margin was 4,585,867,264 bytes (about 4.27 GiB), so the disk preflight passed. The runner repeats the space check before every checkpoint and stops if one estimated checkpoint plus the 5 GiB reserve no longer fits.

## Remaining gate

`execution_authorized` remains false. C must issue an authorization matching the committed config, runner, collector, Phase2 lock, R2 input and R2 structure SHA values before `launch_r2.sh --run` can start. Only this single R2 point is requested.
