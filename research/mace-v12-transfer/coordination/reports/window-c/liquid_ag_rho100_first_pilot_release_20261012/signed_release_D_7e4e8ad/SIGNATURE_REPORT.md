# C execution release review — D rho100 first pilot

Status: **SIGNED, narrowly authorized**. This independent execution authorization leaves the frozen science lock byte-for-byte unchanged with execution_authorized=false.

The release is bound to D branch head 7e4e8ad147928f764b152c7c851171bf1f233493, lock SHA-256 e1adaf12d542184c709d2eb2a5363cbfe9578c168ebc5e7ad260ee39a4faf09b, runner SHA-256 cfe4c104ccfa9c69d5aad0b11287edc36a685a37bf9f67fdd8b792b546150fb0, environment SHA-256 f3f6d27f7e8b387bb5cc7b91048dba7fa26280f9ed54aa1572c73879946e867f, and Ag EAM SHA-256 958f4623e70b51c05e26e9b980f3c85bb1ac9852f3f260c8f75c6181e0c58447.

Independent checks passed: current fetched D branch head matched the reported commit; D manifest, preflight, test report and all entries in D SHA256SUMS matched; preflight status was PASS with zero formal MD steps; eight engineering tests were reported passing; checkpoint/RNG segmented-restoration checks were exact; the fixed lock still limits the pilot to step 0–1000 and 3600 seconds. The runner independently checks the authorization file hash and reruns the live preflight immediately before creating formal output.

Wall-time interpretation: do not begin a new integration step at or after 3600 seconds. A step already in progress may finish and be checkpointed; elapsed time may therefore exceed 3600 seconds by that one in-flight step. No subsequent step may begin.

This authorization permits only C-D-AGL-1250K-rho100-pilot01, at most 1000 added steps, with no automatic continuation. It does not authorize the other five Ag tasks, DFT, MACE training, or promotion of the output to a liquid sample or any TRAIN/DEV/TEST role. Any live hash, resource, identity, or preflight mismatch keeps D on HOLD. After the run, D must create an immutable public-safe archive on its branch, fresh-fetch it, restore it, and verify member and outer SHA-256 values before reporting to C. The potential file and paper PDFs must not be uploaded.

Authorization file SHA-256: 623ba38ec086c59719201eb4b4f99922407f09c74e0098de2dcbf1af95eb1e15

This is a repository-commit attestation from C; no GPG key signature is claimed.
