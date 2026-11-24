# V11 validation artifacts (2026-11-25)

This additive directory records V11 geometry, static/cold MD, mapping-only, and short paired TTM software-link validation. It leaves `v11_generated_data/` and all V5–V10 history unchanged.

**Current status:** the static exploratory screen, 10 Å mapping-only preflight, and short no-source/Ag-Gaussian TTM software-link runs passed their stated checks. The 300 K cold baseline still fails: the Ag/W lattice-temperature split grew from about 83 K to 92.5 K during the 1 ps continuation, and potential energy continued to decrease. The TTM runs do not establish melting or physical validity. The Langevin bath supplied about 32 eV net to the atom system in each 0.08 ps pair run, comparable to the 66.1 eV laser input in the source run.

See [V11_COLD_VALIDATION_REPORT.md](V11_COLD_VALIDATION_REPORT.md) for numerical details and limitations. Reproduce the paired TTM software-link runs with [v11_ttm_pair_commands.md](v11_ttm_pair_commands.md). The archived restarts are diagnostic states only, not approved laser-start configurations. `SHA256SUMS.txt` covers this directory except itself.
