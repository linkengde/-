# A review of the V17 development-validation runner

- Attempt 1 passed the static authorization/geometry preflight but stopped at GPAW import before SCF: GPAW 26.7 reported `Please use gpaw python to run in parallel`. No label was produced; log SHA256 `c8131e80b0fcb3bfbd7f701d07d76142cc8d553fa94ad59819163d3a5df72297`.
- Attempt 2 used `gpaw python` but stopped before SCF because the GPAW CLI parsed the driver options. No label was produced; log SHA256 `0788bbd1b1b1e00e6e3daffaf4cf6cea89a62c1eb96142190afe9e3212bfd8c6`.
- Corrected the guarded runner to insert `--` between the GPAW Python script and its `--workspace-root/...` arguments, and to use the selected venv Python for ASE preflight. The pinned structures/roles/PW-PBE recipe are unchanged.
- Confirmed `gpaw python script.py -- --workspace-root ...` passes application arguments intact; the four-rank communicator and argument-forwarding smoke test is run before the next SCF restart. These are launcher fixes only; no DFT result exists for the failed attempts.
