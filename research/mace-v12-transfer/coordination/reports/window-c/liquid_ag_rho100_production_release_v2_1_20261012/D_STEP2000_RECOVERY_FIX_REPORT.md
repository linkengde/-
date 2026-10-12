# D step-2000 manifest and archive recovery fix

The manifest writer raised the undefined-variable error after the physical step-2000 checkpoint had been written. The writer now receives the immutable run config explicitly, records the checkpoint and checkpoint-binding SHA values in the manifest, and uses the archive backend to generate and verify member checksums.

The continuation path reads the immutable archived step-2000 checkpoint, validates its archive commit, full member manifest, checkpoint and sidecar hashes, source step-1000 checkpoint, potential/environment/parameter locks, QA monitor, and the 1600 K thermal log. The failure marker is corrected only in memory for subsequent control flow: QA remains PENDING, and no physical checkpoint fields are changed. This recovery interpretation is explicitly presented to C for authorization.

Validation: test_production_runner_v2_1.py passed all 20 tests, including manifest/checkpoint SHA, real archived step-2000 restore, and fresh-fetch/restore via a temporary bare Git remote. The live read-only preflight is PREPARED_HOLD; the D branch head and remote matched at the tested release commit. EAM parsing, environment SHA, resource gate, historical Q6 regression, and step-2000 position/momentum/RNG restore all passed. The config still has execution_authorized=false. No force/energy calculation, MD integration, DFT, or training was run.

| Artifact | SHA-256 |
|---|---|
| Runner | 41b60d4bb12420ac2f4cba3b041057358e1423b4a580efcf54564595945def6f |
| Config | e0b8effafecadab3de32da7ae135ee99e2b5ee40ad8bf8a06384febfbb585d31 |
| Criteria | dadf5647cd1ffb68fd421e8f2742117c604b324e20443c0e2551e7208ae21c69 |
| Step-2000 checkpoint binding | b314ca230c412c5852e289e755dcda54ceefd0746d8dd2eff6bf9d1bcb56a285 |
| Step-2000 checkpoint | 3ae061355544c7c1871605d11cc532c58d975fd0d77611dfa6fab26b2f238381 |
| Environment | e71e04b8edc989001a702ffb308b1b897d6b521d0ebfdcb41b43592ae2bff1fb |

D code release commit: 2fb55ac472b2175f39bbe4f18edcd426077e413d. The C authorization request is in C_STEP2000_CONTINUATION_AUTH_REQUEST.json. Do not start continuation until C signs an authorization bound to the current D remote HEAD and all hashes above.
