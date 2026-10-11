# Ag64 single-point preparation

This is a future execution recipe for a C-reviewed liquid Ag candidate. It does not create a structure and is not an authorization to run GPAW.

## Entry conditions

1. Receive a frozen candidate structure and source record. Record the source/trajectory lineage, sampler and its reviewed status, temperature, density and density reference, seed, parent-family ID, and candidate selection rationale.
2. Record the task ID, input file SHA-256, atom count (64), element order/IDs, cell and PBC. The structure must be an existing reviewed input; do not synthesize or repair coordinates in the DFT run script.
3. Record a C-reviewed execution release and its exact configuration SHA. Check that the task is included in the release and that the F host is released.
4. Confirm the chosen storage path is durable, has write permission, and meets the checkpoint and archive budget. The current `/workspace` overlay does not meet this condition by itself.
5. Freeze the method record: GPAW 26.7.0, ASE 3.29.0, PBE, gpaw-data 1.2.1, PW cutoff 500 eV, Fermi-Dirac electronic smearing 0.10 eV, and energy/density/eigenstate convergence thresholds `1e-5`. Keep the ionic sampling temperature separate from electronic smearing.

## Same-geometry k-point gate

For liquid Ag, register the exact same input hash for a Γ calculation and a 2×2×2 calculation. This is a numerical comparison, not two candidate structures. Compare total energy per atom and per-atom force vectors. The handoff's pre-registered limits are ΔE ≤ 2 meV/atom, force-difference RMS ≤ 0.01 eV/Å, and maximum force difference ≤ 0.03 eV/Å. If any limit fails, register and review a 3×3×3 or denser comparison before selecting a production mesh. Preserve each mesh result as a control record; do not silently use the failed coarse mesh.

## Fixed-ion single point

Once the structure, exact method, mesh and release are frozen, run one static-ion calculation. Use the ASE/GPAW route already present in F's environment. Keep the input's coordinates, atom order, cell and PBC byte-for-byte/array-for-array unchanged. Use the handoff's mixer and iteration cap only if they are included in the reviewed per-task manifest. Do not relax ions or cell as part of this single-point step.

The launch form for a later, separately released task is:

```text
mpiexec -n 4 gpaw python calculation.py --manifest approved_manifest.json
```

Do not run this command during preparation. A template is not a release.

## Required checks after a future run

- SCF reports `converged=true` and does not exceed the approved iteration/time budget.
- Native GPAW energy, force-consistent free energy, all 64 three-component forces, and the full stress tensor are finite.
- Output structure has exactly the input elements, order, IDs, coordinates, cell and PBC.
- The result verifier confirms the input hash, method/PAW hashes, mesh, convergence, label-key convention and checksums.
- If the label is intended for a thermally sampled training candidate, pair `free_energy_eV_cell` with the matching forces as directed by C's frozen data-role manifest. Keep native extrapolated energy as a separately named diagnostic; never merge the two energy conventions in one target column.

See `AG64_INPUT_MANIFEST.template.json`, `RESULT_FORMAT.template.json`, and `ARCHIVE_LAYOUT.md`.
