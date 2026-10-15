# Stage C-D5 — M77 geometry and admission audit

This package contains three geometry-only 77-atom Ag/Ti₃SiC₂ interface candidates and an admission audit. No DFT, LAMMPS `run 0`, NVE/NVT, TTM heating, fitting, or potential edit was performed.

## Files

- `MINIMAL_MODEL_COORDINATES/Ag_Ti3SiC2_M77_gap_2p16A.extxyz`
- `MINIMAL_MODEL_COORDINATES/Ag_Ti3SiC2_M77_gap_2p40A.extxyz`
- `MINIMAL_MODEL_COORDINATES/Ag_Ti3SiC2_M77_gap_2p64A.extxyz`
- `ATOM_ORDER.csv`: stable 1–77 atom ID, element, region, layer, and site mapping.
- `MODEL_MANIFEST.json`: geometry parameters, source hashes, structure hashes, and geometry metrics.
- `POTENTIAL_AND_RUNTIME_GATE.json`: actual LAMMPS/runtime and potential identity gate results.
- `STATIC_GEOMETRY_METRICS.csv`: distances and the explicit unavailable `run 0` force gate.
- `SHORTEST_20_DISTANCES.csv`: shortest 20 periodic pairs for each candidate.
- `generate_minimal77.py`: deterministic geometry-only generator, with source file SHA recorded in the manifest.
- `STATIC_ADMISSION_REPORT.md`, `300K_STABILITY_REPORT.md`, and `D4J_MODEL_SELECTION.md`.

The substrates use the seven-plane sequence Si–Ti–C–Ti–C–Ti–Si. The 2.40 Å gap is a geometric reference, not a relaxed equilibrium spacing. The coordinates have a 0° aligned area-7 cell and do not retain the parent 12.51983° twist.

`D4J_MODEL_READY=NO`: the requested LAMMPS force calculation and 300 K NVE/NVT checks remain blocked because no LAMMPS runtime is available here and the only found EAM+Tersoff+Morse combination is explicitly an unvalidated calibration-v0 candidate.

Recreate in a repository checkout with ASE 3.29.0:

```bash
/workspace/.venvs/mace/bin/python generate_minimal77.py
```

