# Ag64 rho100 pilot run result

- Task: `C-D-AGL-1250K-rho100-pilot01`
- C authorization commit: `83f67ff2ef3799fb3893038e88df5de4f393cab6`
- Authorization SHA-256: `623ba38ec086c59719201eb4b4f99922407f09c74e0098de2dcbf1af95eb1e15`
- Scientific lock SHA-256: `e1adaf12d542184c709d2eb2a5363cbfe9578c168ebc5e7ad260ee39a4faf09b`
- Runner SHA-256: `cfe4c104ccfa9c69d5aad0b11287edc36a685a37bf9f67fdd8b792b546150fb0`
- Environment SHA-256: `f3f6d27f7e8b387bb5cc7b91048dba7fa26280f9ed54aa1572c73879946e867f`
- EAM SHA-256 (file not included in this archive): `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`
- System: 64 Ag atoms, fixed rho100 cell, periodic in all directions
- Integrator: ASE Langevin NVT, target 1600 K, timestep 0.5 fs, friction 0.01 fs^-1

## Completion

- Actual steps: **1000** (`0 → 1000`)
- Simulation time: **0.500 ps**
- Wall time: **197.483 s** (3600 s maximum)
- Stop status: **COMPLETED_MAXIMUM_STEP**; normal stop at the authorized step limit: **YES**
- Failure: `None`; automatic continuation: **NO**
- Final checkpoint SHA-256: `c820437a2512c4b3d7dc9e60bb60ac3ced6588aa1f2a629d58cfdabb24f34046`
- Trajectory frames: 11 at steps 0, 100, …, 1000; checkpoint contains positions, velocities/momenta and all 3 PCG64 states.
- Final minimum periodic pair distance: **2.424554 Å**; all 64 Ag atom IDs are unique and all checkpoint coordinates and velocities are finite.

## Logged thermal and energy summary

The logs contain 51 samples at 20-step intervals. Temperature statistics below are descriptive of this short hot-start pilot, not an equilibrium assessment.

- Target temperature: 1600 K
- Logged temperature mean: **2258.446 K**
- Logged temperature sample standard deviation: **637.012 K**
- Logged temperature range: **1597.558–4086.011 K**
- Total energy at start/end: **-100.478344 / -144.831780 eV**
- Endpoint total-energy change over 0.5 ps: **-44.353437 eV** (**-0.693022 eV/atom**)

This is an NVT Langevin run, so total energy is not conserved. The logged mean and peak temperature substantially exceed the 1600 K target during this short trajectory; the energy endpoint change is also reported without treating it as a conserved-energy drift. These observations require C's scientific review. This pilot does not establish an equilibrated liquid, and it is not assigned a TRAIN/DEV/TEST role.

## Integrity and scope

The runner's final checkpoint SHA agrees with `checkpoint.sha256`; its stored hashes for `trajectory.traj`, `thermal.csv`, and `energy.csv` match the current files. The immutable archive excludes `Ag.eam.alloy`, PDFs, credentials, and unrelated tasks. No second trajectory, DFT, MACE training, or continuation was run.
