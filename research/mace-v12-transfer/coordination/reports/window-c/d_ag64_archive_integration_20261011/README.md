# Stage D: production GitHub archive integration

This stage preserves the frozen step-1000 source archive and the old core adapter. `segment_core.py` is a byte-for-byte copy of the reviewed ASE 3.29.0 integrator. `production_segment_adapter.py` is the new bounded entry point; `github_archive_backend.py` implements append-only GitHub `upload()`, fresh-clone `verify()`, and safe `restore()`.

The only production scope is task `C-D-AGL-1250K-rho097-seed01`, global steps 1000 through 2000 inclusive (1000 added steps). The source archive, potential SHA, old integrator identity and all physical parameters are pinned. Every 100-step bundle must be uploaded, freshly cloned, verified and restored before the next integration step. A push/readback failure stops before another step. The original step-1000 checkpoint is read-only.

The integration test uses a copy of the archived step-1000 atom/RNG state and a deterministic synthetic state-only transformation. It does not construct a calculator, evaluate energy/forces, call ASE integration, or modify the source checkpoint. The test archive is clearly marked synthetic and is not a scientific result.

No step beyond 2000, 120000-step continuation, potential change, checkpoint overwrite, temperature change, RNG reset or velocity reinitialization is permitted. The executable release JSON is generated only after the real GitHub upload/fresh-clone/restore test passes and is hash-bound to the code and exact D branch head.
