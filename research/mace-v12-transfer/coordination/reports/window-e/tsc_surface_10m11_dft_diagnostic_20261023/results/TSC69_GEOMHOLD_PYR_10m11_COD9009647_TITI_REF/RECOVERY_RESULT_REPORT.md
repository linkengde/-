# TSC69 fixed-geometry DFT recovery result

## Outcome

`STATUS=SUCCESS`; GPAW reported SCF convergence after 9 resumed iterations, for 69 total SCF iterations. The calculation resumed from the read-verified step-60 checkpoint and did not change the structure or scientific settings. The checkpoint used for resumption had SHA-256 `5cd9f365ac2b80e888e3f357cb93a12798f59c263cd6b4b823017c16a60d6c2a` and was removed only after the converged result had been written.

## Locked calculation

- GPAW 26.7.0; ASE 3.29.0; PBE/PAW; plane-wave cutoff 500 eV.
- Standard Monkhorst–Pack `6×6×1`; Fermi–Dirac width 0.10 eV; non-spin-polarized.
- Davidson `niter=2`; MSR1 mixer and SCF thresholds from the bound diagnostic configuration; 4 MPI ranks, 1 OMP thread.
- Fixed geometry; no relaxation or other structures were run.
- Configuration SHA-256: `15ddc8f596281298656a5a21df67b2915f0aebd5ba3bb18316fd2935ee8729c2`.
- Input SHA-256: `4b84cfc8f76827c7f81d9401f60e0c161355d61adc57cfe9ba867a78c780f473`.
- Geometry SHA-256: `fd3351b1ba798ea9e0567013740fba8651d801fa8b62199b89da2c4d5ef37669`.
- Recovery-v2 runner SHA-256: `b5671637e47a4a179cf17314212fc57a9909a4bc0a65ee5c021fa6d488eaf299`.

## DFT labels

- Total energy: **−532.2613287689045 eV**.
- Force-consistent free energy: **−532.7984214555453 eV**.
- Forces: 69 atomic vectors, in eV/Å.
- Stress: available as a full `3×3` tensor in eV/Å³ in `RESULT.json` and the labelled extxyz.
- Maximum force magnitude: **2.6255850886 eV/Å**, on source atom ID 63 (Ti).

The large maximum force is retained as a property of this frozen geometry. This is a converged fixed-geometry electronic single point, not a relaxed structure or evidence of mechanical stability.

## Integrity and role

The labelled extxyz was read back and compared with the input: atom count, element order, source atom IDs, coordinates, cell, and PBC match exactly. Result shapes and finite values were checked. The manifest covers all 38 artifacts (excluding the manifest itself), and every entry passes `sha256sum -c`. No checkpoint remains after successful result persistence, and no GPAW process remains active.

`DATA_ROLE=UNASSIGNED_SAME_SOURCE_FACET_DIAGNOSTIC`; this result has not been assigned to TRAIN, DEV, or TEST and must not be treated as an independent validation family.

## Key output hashes

- `RESULT.json`: `d8b764cf6757f4e1476624445a64b0deb5e63ecec412b52a57c62441d043a13b`
- `labelled.extxyz`: `2911d08e576c1967cf8c81cdb89510224c8967266ab27aeea599c9cce0f64a44`
