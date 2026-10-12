# Runner/preflight checks

- Python syntax compilation: PASS.
- Read-only source/environment/resource preflight: `PREFLIGHT_PASS`; `execution_authorized=false`; no DFT started.
- Missing-authorization gate test: PASS; `--run` exited 2 with `C execution authorization is missing; HOLD` before GPAW imports or SCF construction.
- Requested calculation: not run in this authorization package.
- Runner SHA-256: `36f35f1657a8fb9ae80c5b6dc4bf42869b694ba891f496a20739626826613f5d`.
- Config SHA-256: `03dcd34d74249c360fc16414c48c3d5fe53d175bc855248a0e3403e6ee3f17b8`.
