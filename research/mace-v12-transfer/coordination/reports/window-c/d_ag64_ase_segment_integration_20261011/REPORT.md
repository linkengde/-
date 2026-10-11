# ASE 3.29.0 step-1000→2000 adapter integration review

**Decision:** `CODE_REVIEW = PASS_FOR_SYNTHETIC_ADAPTER`; `NEW_MD_RELEASE = NO`; `D_CONTINUE = HOLD`; other five tasks, DFT, and MACE training remain HOLD. The code SHA is frozen for review, but production GitHub archiving is not integrated or verified. No execution instruction was sent to D.

## Source checkpoint verified read-only

The D archive from commit `145ee8050cf95a3a1988e868c6b935e10538f723` has archive SHA-256 `b62a818e58a2cf423e4cb1b9ed0468ae73ad578a75c4b2cb9fb3f28baf75e7b1` and external manifest SHA-256 `ab9f0dffa8ef145eac7792dfe1e14e07c1354df1c2152dafec97bad1dc3a57b1`. All 18 members match the external manifest. The step-1000 state file, array file, and source manifest match their inner hashes.

The adapter reconstructed 64 Ag atoms, positions, momenta, masses, cell, PBC, atom IDs, global step 1000, and the packing/velocity/thermostat PCG64 states. Positions, momenta, masses, cell, PBC, IDs, and reconstructed velocities compared bit-for-bit with the archived arrays. No calculator was attached and there were zero energy, force, or integration calls. See `SOURCE_RESTORE_PROBE.json`.

The archived `velocities_A_per_fs` array label does not match its stored unit convention: its values equal ASE's raw `Atoms.get_velocities()` and `momenta/masses`. The adapter treats momenta as authoritative, performs no velocity conversion, and requires exact reconstructed velocity equality. New bundles use the explicit name `velocities_ASE_internal_units`.

## Hash scope and restart semantics

The old `scientific_parameters_sha256` recomputes exactly to `8f792d6206c1b7bae7b900e7ea81445d6cd3e837de2b5979c88bb977650ab7b2`. Its `fixed_identity` includes potential identity, non-resource sampling fields, task row and source, package versions, old runner/parser code hashes, work root, ASE Langevin version, host/Python/thread metadata, and thermostat description. It is a historical provenance digest, not a pure physics-only identity or the new code hash. The original runner SHA `6b54d094c26bb9d4a5e4eb5997a31955eb8829a2975042e23fe10527cb943583` was independently recomputed and matches its manifest record. The new adapter SHA is kept separately.

ASE 3.29.0 `Langevin.step()` draws `xi` then `eta` from the supplied Generator, updates positions, recomputes forces, and updates momenta. `set_temperature()` recomputes coefficients without consuming RNG. The original uses `FixCom` with `fixcm=False`; this adapter restores raw momenta before attaching the same constraint and does not project them during load. Time step remains `0.5 fs`, friction `0.01 fs^-1`, and the global schedule selects 1600 K for steps 1000–1999. The original source continues to set its default continuation endpoint at 120000; the new adapter does not invoke it.

The original thermal log's `sampler_energy_eV` is potential energy (`get_potential_energy`), not total energy. Kinetic energy is a separate quantity, and original pressure includes the kinetic correction. New segment logs keep potential and kinetic terms separate and preserve the same pressure expression.

## Integration and synthetic test results

`ase_segment_adapter.py` subclasses the real ASE 3.29.0 integrator. Its overridden `step()` checks the hard global stop before ASE can calculate forces or draw random numbers, performs one ASE integration step, validates finite arrays, and increments the global step inside that same method. `run()` and `irun()` are disabled. The output loop writes 20-step thermal rows, 100-step trajectory frames, and immutable 100-step state bundles in that order; an archive/readback failure stops before the next step.

All 14 synthetic tests passed. Continuous and interrupted/resumed runs matched bit-for-bit across all 1000 per-step position, momentum, and PCG64 states. A step-2001 attempt was refused before force or RNG use. Tests injected timeout, nonfinite force, wrong task/owner/boundary/code hash, unsigned lock, remote gate failure, archive failure, and existing-output collisions. Local synthetic tar fresh-extract and internal SHA restore passed. The synthetic output contained 50 thermal rows and 11 trajectory frames through step 2000. Full results are in `SYNTHETIC_TEST_REPORT.json`.

Only a marked virtual harmonic calculator was used for the integration tests. The D checkpoint itself was only restored and inspected; no EAM potential was read or evaluated. No real MD, DFT, or training ran.

## Frozen review hashes and remaining blockers

- `ase_segment_adapter.py`: `e0bba3702c54bbdcb39f5f7e5641e1e1969a2ea5af8faceee8a0d50316f440dd`
- `test_ase_segment_adapter.py`: `9185ed5f4e5b4410fb2568e51d2789d096941e53a3d8696cc05bf82673e60de3`
- Candidate lock: unsigned; its SHA is listed in `SHA256SUMS.txt` and the final commit pins all package files.

Two hard blockers remain: (1) no real GitHub branch push/fresh-fetch/restore test for a new segment archive from D; (2) the CLI deliberately refuses execution because that archive backend is not implemented. D's actual fresh-fetched archive paths and live branch head also need a new read-only check. The current command draft is explicitly blocked and must not be sent as an instruction.

The candidate preserves the source physical values: 64 Ag, 9.032733418760 g/cm³, 10.826824658074 Å cubic cell, seed 125009701, 1600 K hot segment, 0.5 fs, and 0.01 fs⁻¹. The bounded segment is only global step 1000→2000; the existing 120000-step path is never used. See `CANDIDATE_RELEASE_LOCK_UNSIGNED.json` and `LAUNCH_COMMAND_DRAFT.md`.
