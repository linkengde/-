# Ti₃SiC₂ backbone and Ag/Ti₃SiC₂ interface DFT input package

This package prepares fixed-geometry DFT inputs from the E7 audit. It does **not** contain DFT labels, run DFT, modify TSCBOP, or start fitting.

## Prepared inputs

- Six Ti₃SiC₂ bulk points: ±2% volume, ±0.5% shear, and a local Si z displacement of ±0.02 Å. These are copied byte-for-byte from the E7 geometry-only package. All belong to `TSC_KISI_COD9009647_BASELINE_DERIVED_FAMILY`.
- Nine Ag/Ti₃SiC₂ Si-termination interface geometries at plane gaps 1.8–3.4 Å in 0.2 Å increments. Each comes from one frozen 217-atom parent and differs only by rigid translation of the Ag slab along z. All belong to `AG_TSC_SI_TERM_TWIST12P51983_217_PARENT_FAMILY`.
- Per-structure JSON DFT configurations, a shared authorization-gated fixed-single-point runner, input manifest, family ledger, and family-level split proposal.

Every item is `TRAIN_CANDIDATE`/`GEOMETRY_ONLY`; its DFT labels do not exist. The runner requires a separate authorization file binding task, configuration SHA-256, geometry SHA-256 and runner SHA-256. No such authorization is included. Running it without an exact authorization exits before GPAW is imported.

## Method and limits

The six bulk points use the E6/B-matched GPAW PBE/PAW PW-500 eV, Fermi–Dirac 0.10 eV, Pulay β=0.05 (`nmaxold=8`, `weight=100`), 1e-5 convergence criteria, 160 SCF iterations, four MPI ranks/one OpenMP thread, and the existing 10×10×4 mesh. Their configs bind the audited Ti/Si/C PAW setup hashes. They request native energy, free energy, forces, and full stress.

The interface scan uses the previously proposed Gamma-centered 2×2×1 mesh with the same functional, cutoff, smearing, mixer and convergence thresholds. This mesh has **not** been shown converged for this 217-atom cell; the inputs are not cleared for DFT execution until a method/resource review. The slab's 3D stress includes vacuum dilution and is diagnostic only, not a bulk stress target. The parent is unrelaxed and starts at a 2.0 Å plane gap; the scan does not claim that any point is an optimized equilibrium interface.

The 217-atom coordinates are from the archived geometry-only commit `f1a080e14cd5ffd9739dcd2a3df2d188c853a8a8`, file SHA-256 `945c5e214c1ceaeba4125ff785f4092e14858c75a985061a7841da3e5ce6a9a2`. Its metadata records a 394-atom parent SHA, but the full parent file and generation script are not part of this package, so deeper ancestry reproduction remains incomplete.

## Split decision

Keep each parent family indivisible. The E7 baseline-derived family and the interface scan are proposed TRAIN families only. Existing baseline/Ti/C reference labels remain numerical-control/diagnostic and are not silently transferred into TRAIN. No independent DEV source was found in E7; TEST remains evaluator-controlled and was not read or enumerated. Therefore this is an input package, not a fit-ready dataset, and it cannot support DEV-based model selection or an independent generalization claim.

The interface data can inform a future unified Ag–Ti–Si–C model. Ag-containing labels cannot constrain three-element Ti–Si–C TSCBOP parameters. The two objectives must remain separate.

## Reproduce the input construction

```bash
/workspace/.venvs/gpaw-mpi/bin/python build_input_package.py
```

This copies the six E7 structures, rigidly shifts Ag in the nine interface geometries, and writes the configs/manifests. It performs geometry I/O only; it does not instantiate a calculator. Validate all bytes against `SHA256SUMS.txt` after regeneration.

`run_fixed_singlepoint.py CONFIG.json AUTHORIZATION.json` is a draft execution entry point. Its expected failure without an external signed authorization is intentional. It is not a DFT release.
