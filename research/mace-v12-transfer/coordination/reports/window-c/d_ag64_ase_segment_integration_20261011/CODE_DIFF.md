# Adapter integration diff

- Added `ase_segment_adapter.py`; left the original sampler and prior proposal untouched.
- Added a bounded `Langevin` subclass. Its actual `step()` checks global step `< 2000` before temperature changes, RNG draws, force calculation, or coordinate updates. It increments the global step internally only after finite state validation. `run()` and `irun()` reject calls, so ASE's default endpoint cannot bypass the bound.
- Restart reconstruction uses the archived arrays and full PCG64 states; no Maxwell-Boltzmann redraw, reseed, velocity initialization, or position generation exists in the adapter.
- Added immutable per-boundary bundle and tar creation. Output includes a state NPZ, full RNG state, segment trajectory frames, thermal log, README, MANIFEST, and SHA256SUMS. Each member and tar archive is capped at 50 MiB. Existing names are refused and publication uses no-replace directory rename / hard-link archive creation.
- The callback order is integrated step and finite-state check, thermal row every 20 steps, trajectory frame every 100, then checkpoint/archive at each 100-step boundary. A failed remote callback raises before another ASE step; the source step-1000 files and last local immutable bundle remain.
- Actual GitHub remote orchestration is intentionally not implemented in the CLI. It is a hard release gate and the CLI returns HOLD rather than running local-only scientific work.
