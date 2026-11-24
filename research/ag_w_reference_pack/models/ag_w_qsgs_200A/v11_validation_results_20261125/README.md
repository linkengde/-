# V11 validation artifacts (2026-11-25)

This directory records a local Dai EFS static and 300 K cold-baseline run on the committed V11 bicontinuous Gyroid geometry. It is kept separate from `v11_generated_data/` and all V5–V10 records.

**Result:** exploratory static screen passed, cold baseline gate failed because the Ag/W lattice-temperature difference remained about 83 K and potential energy was still decreasing. A separate 1 ps NVT continuation did not resolve the phase-temperature split. No real-geometry mapping or TTM/source run was performed.

See [V11_COLD_VALIDATION_REPORT.md](V11_COLD_VALIDATION_REPORT.md) for setup and interpretation. The 1 ps continuation script is `in.v11_cold_followup`; it resumes the archived step-8200 restart. The archived restarts are diagnostic only; do not treat them as an accepted laser-start state. `SHA256SUMS.txt` covers files in this directory except itself.

To decompress the archived cold snapshot before loading it in LAMMPS:

```bash
gzip -dc v11_cold.restart.gz > v11_cold.restart
```
