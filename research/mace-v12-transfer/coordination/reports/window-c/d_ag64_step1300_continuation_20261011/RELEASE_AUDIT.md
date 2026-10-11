# Release audit

## Decision

`D_CONTINUE = AUTHORIZED_FOR_THIS_SEGMENT_ONLY` in the signed configuration. C did not start or dispatch MD. D's live preflight is still required immediately before execution; if the remote D branch head, checkpoint, potential, ASE version, or output-root condition differs from the signed pins, stop and report.

## Verified evidence

- D branch head: `9aea23179dd73da749d67f7be0b3a3b86e5df085`.
- D archive receipt SHA-256: `ede5dfee81846911e289540391432c5dfe324aac159cd89aa6dc08cafb93335f`.
- Archive SHA-256: `e75322b9bdd4b6bc2d485225a76c889536219ac22b6fe299858f96c1411839ee`.
- Checkpoint state SHA-256: `3b3950729735765dc5c4d205d16a3654b4e9f8ba193c58bf8421685f36d66b9e`.
- State arrays SHA-256: `6eb03c97a33f6c2f3b5d9521c98435a77e5ee9228286090602733c3d22a65a2a`.
- Potential SHA-256: `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`.
- Existing archive backend SHA-256: `8fd001a9b1706fc8832e90a77c22783a22a77eab06f50a66bd1d1eeb97c00346`.
- Prior step-1000 C release SHA-256: `026be260ca65ab55b60c96a407c6181649e5f245bcf2d138993e246d2aa56d9b`.
- Offline test runtime: Python 3.12.14, ASE 3.29.0, NumPy 2.5.3.

The C test fresh-cloned D's branch, verified its pinned head and archive receipt, checked all archive and bundle digests, reconstructed the 64-atom state and RNG streams, and confirmed the integrator refuses at the 2000 boundary before accessing a calculator. No forces or integration steps were evaluated.

## Preserved limits

The signed range is only 1300 to 2000, at most 700 added steps. The original potential, positions, momenta, velocities, RNG streams, cell, PBC, atom ordering, temperature path, timestep, thermostat friction, seed, atom count, and density are pinned. The 600-second wall limit is unchanged. A timeout or archive failure stops work without automatic continuation. D must use a new non-overwriting output directory and keep every checkpoint archive recoverable.

The raw source `execution_status` marker is retained and disclosed in the signed config. Its provenance binds the prior C-signed release and declares `synthetic=false`; the user has reported that the D interruption audit is complete. The discrepancy is recorded without altering the prior archive.
