# Release audit: step 1600 to 2000

- Task: `C-D-AGL-1250K-rho097-seed01` (D only).
- Authorized global interval: 1600 through 2000; maximum 400 added steps; step 2001 and automatic continuation to 120000 are forbidden.
- Checkpoint source: D branch `cloud-d-b1-dft`, head `9ff9cee57714d00128607828d6e5afd0eaf6e946`.
- Archive SHA-256: `55f646d7c2b3069a6c54ae8e9bf54c21b65a2fe0d0239067699a14239275b1c0`.
- Archive receipt SHA-256: `4cbd28618d82e3b435b7e1bdf1ffa34ea638fae474e882aa838df71f0807a87f`.
- `state.json` SHA-256 / supplied checkpoint SHA: `3688c9e9b0b3587e9c7ae275be32e7366eb9236485404cd7cb236993118fbba6`.
- `state_arrays.npz` SHA-256: `35c8256ea58b9e7d44033e5cdcfb9ef3b515de822bafdb919f03ad3cf7829701`.
- Restored identity: 64 Ag atoms, global step 1600, PBC and cell accepted by the pinned reader, RNG streams `packing`, `velocity`, `thermostat` preserved as PCG64 state.
- Potential SHA-256: `958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447`.
- New adapter SHA-256: `0ba600e063b11af179343d1cd2dd92b891736aa75d68f5c2515a97b8b8bff7b4`.
- New continuation core SHA-256: `53b51c7ec048e7fd15ce101431a1ddfceff59efdc2b68ed2948e2d399ff7e082`.
- Archive backend SHA-256: `8fd001a9b1706fc8832e90a77c22783a22a77eab06f50a66bd1d1eeb97c00346`.
- Signed config SHA-256: `40499271b9d808fe7e34ec737d96599df5c742ffb5a409d0d456c97d16aab34c`.

Verification performed: fetched D branch head pinned above into a fresh temporary clone, verified archive receipt and outer archive digest, validated member and manifest SHA-256 values, restored state arrays without attaching a calculator, checked global step and RNG stream identity, and validated the signed release lock against the adapter/core/backend hashes. No integration step or force calculation was run.

The archived step-1600 provenance identifies its source as the previously authorized 1300-to-2000 continuation, with `synthetic=false`. That original provenance is retained. The new release does not alter the archive or prior release package.

The inherited wall-time cap is 600 seconds. The previous 1300-to-1600 segment used about 605 seconds for 300 steps, so the 1600-to-2000 target may be interrupted by the wall-time gate before reaching step 2000. Any such partial result must remain HOLD and require a new release based on the latest verified checkpoint.

`execution_authorized=true` means the signed config is prepared for D's single bound continuation; C did not start or dispatch it.
