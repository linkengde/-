# Execution readiness — Ag–Ti/Si/C minimum DFT preparation

```text
READ_ONLY_AUDIT_COMPLETE_FOR_ACCESSIBLE_SIX=YES
REQUESTED_E_HANDOFF_RECEIVED=NO
135_MATCHED_COMPONENTS_VERIFIED=NO
EXISTING_RECORDS_VERIFIED=6
EXISTING_PARENT_FAMILIES_VERIFIED=3
NEW_STRUCTURE_COUNT_FROZEN=NO
QE_RUNTIME_READY=NO
QE_PSEUDOPOTENTIALS_VERIFIED=NO
PEAK_MEMORY_BENCHMARKED=NO
DFT_EXECUTION_READY=NO
START_DFT=NO
```

The six accessible historical records pass SHA, geometry/ID, complete-force, archived SCF and PAW-fingerprint checks. They can retain their existing diagnostic use without repetition. The E-reported 135-atom matched components and its at-least-ten-new-SCF plan have not been received in C and cannot yet be verified. The sender's checksums and counts are recorded as sender assertions, not C passes.

The current memory snapshot has about 31.64 GiB headroom; the previous 5.62 GiB observation is stale for this task. Resources are nevertheless **not ready**: QE `pw.x`, the UPF set, actual new inputs, measured memory/scratch requirements, and validated QE restart storage are missing/unverified.

## Minimal next actions

1. Transfer the existing 634,192-byte handoff tar.gz to this C conversation, or publish its already prepared files at a reachable Git commit/path. Expected tar SHA: `4e35ce85e6a47295e07c6e4e688925bc343fee7cdd966d39e3619b9d00d7026f`. C's local `/workspace/-` is not E's local filesystem. Do not regenerate existing structures or repeat SCF to solve a transfer problem.
2. Verify the tar bytes and internal 208-file manifest; read `DFT_REQUEST_TO_C.md` (if supplied) and `C_MINIMUM_DFT_TASKS.md`, then audit matched 135 full/phase results and deduplicate the V19/old-six overlap.
3. Freeze the actual minimum task matrix and identify whether E's ten-point lower bound still applies; document contact-channel coverage and independent validation provenance. The conditional one-parent gap plan is 2–3 new full geometries and 4–7 new SCF inputs depending on a reusable parent triplet; it is not the complete three-channel fitting requirement.
4. Review the exact QE software/UPF/input method and arrange a bounded resource/convergence pilot only after a separate formal calculation instruction. Do not replace QE with GPAW or alter precision without review.

No TSCBOP, Ag EAM, Morse parameter, existing data role, historical Stage67 result or sealed record was modified. No new scientific structure was generated. This delivery is a partial evidence-backed review and remains HOLD pending the existing handoff, not an execution release or a fitted-potential result.
