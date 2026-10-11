# D Segment Archive Backend Synthetic Test

This directory contains inert synthetic bytes for testing archive upload and restoration only. The mock checkpoint and mock trajectory are not real scientific artifacts.

Test route: local package and deterministic tar archive -> normal push to the D-owned GitHub branch -> independent fresh clone -> extract into a new temporary directory -> compare SHA-256 for every package member.

No step-1000 checkpoint, real potential, real trajectory, or sampler state is included or read by this test. No MD, DFT, or training is run.
