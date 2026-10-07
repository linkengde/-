# V17 parent hash and energy-convention audit

This is a read-only provenance audit of the three assigned V14 training parents
(AgC frame 29, AgSi frame 30, AgTi frame 28) and the eight already-published V17
paired training-diagnostic archives. It does not access V17 development-validation
outputs or the six withheld-test labels, and it launches no DFT or model work.

## Reproduce

From the repository root, using the existing GPAW/ASE environment:

```bash
/workspace/.venvs/gpaw-mpi/bin/python \
  research/mace-v12-transfer/coordination/reports/window-b/v17_energy_convention_and_parent_hash_audit/audit.py \
  --repo-root /workspace/-
```

The script verifies source and archive hashes, original GPAW settings and log
energies, the exact REF energy/force copies, the V17 canonical parent hashes,
and each of the eight child input displacements. It writes `audit.json` and
`parent_evidence.csv` beside itself. The historical V15 hash algorithm is not
reconstructed because its generating script and serialization contract are not
present in the committed source tree.

See `report.md` for evidence interpretation and the V17 integration recommendation.
