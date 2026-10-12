# Ag64 rho100 pilot release package

This C-only package freezes the scientific identity and parameters for one D-owned Ag64 rho100 exploratory segment. It does **not** authorize execution yet. `EXECUTION_AUTHORIZATION.json` remains `false` until D publishes the missing live engineering/preflight evidence and C binds the exact hashes in a new authorization revision.

The source PDFs and the Ag EAM parameter file are not included. Only their citations and SHA-256 values are recorded. The 1000-step segment is a 1600 K warm-start run, not 1250 K equilibration or proof of a liquid state.

Files:

- `LIQUID_AG_PILOT_LOCK.json`: frozen density, cell, potential identity, initial packing, RNG streams, dynamics and stop/restart contract.
- `EXECUTION_AUTHORIZATION.json`: current explicit HOLD record.
- `REPORT.md`: evidence review and remaining release gates.
- `SHA256SUMS.txt`: package integrity list.
