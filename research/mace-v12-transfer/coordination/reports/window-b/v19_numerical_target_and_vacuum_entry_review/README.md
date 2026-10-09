# Numerical target and disabled vacuum/dipole entry review

Six compact 5x5/6x6 archives at Fermi widths 0.05/0.10/0.20 eV independently pass their registered verifiers and full SHA256 inventories. Eight comparisons are recomputed without writing production files. At width0.20, 5/6 vector difference0.0099259526 eV/A passes the unchanged0.01 budget, energy0.272793747 meV/atom and pair0.007120371 eV/A also pass. This is one local comparison, pending A6/8; SCF convergence and broader occupation smoothing do not establish zero-temperature accuracy. Width0.05 5/6 vector0.019478066 fails. At mesh6, width0.10→0.20 changes vector forces0.030257645 and pair0.047656030 eV/A: occupation target sensitivity is material.

## Energy and force convention

Installed GPAW26.7 source is authoritative evidence for this version, with source hashes in report.json. `new/ase_interface.py:get_potential_energy` returns `free_energy` when force_consistent=True and `energy` otherwise; `new/calculation.py` stores total_free and total_extrapolated separately. ASE3.29 atoms.py documents free energy as force-consistent. `occupations.py:FermiDiracCalculator` uses extrapolate_factor=-0.5; `new/ibzwfs.py` applies that factor to the entropy energy. The legacy equivalent in `old/hamiltonian.py` is E_extrapolated=F+(-0.5)*(-TS). The extrapolation is a finite-width estimate, not an exact zero-temperature calculation. Returning a different energy key does not change computed forces. Native energy paired with finite-width forces must retain its convention and must not be described as an exactly force-consistent target.

Common-target proposal: freeze geometry/PAW/PBE/cutoff/mesh/width/boundary/dipole and energy key as a target identifier before generating references. For a force-consistent finite-width fit use F and its forces at one explicitly chosen width; retain native extrapolated energy separately for diagnostic comparison. For a zero-temperature intended fit first demonstrate width/mesh convergence and energy-gradient consistency; finite-width pass alone cannot select that target. Width and mesh remain undecided pending A6/8. Do not mix sigma0.20 slab energies with Gamma/sigma0.10 historical energies under one undifferentiated target. Historical Gamma/sigma0.10 native rows remain unchanged and identified as legacy; they can be used as a separate declared legacy task or require matched relabeling before a common-target integration.

Relabel requirements: all retained rows with incompatible width/energy convention; periodic slab rows with inadequate Gamma sampling; unresolved legacy rows whose original method is not evidenced; and frame0 with the historical1x4x2/Gamma discrepancy. Relabel exact original geometries, archive F/native/forces and complete method; rebuild a candidate only in a separate authorized integration. Isolated finite clusters do not automatically require the slab mesh. Existing numerical pilots are not training/development labels. No relabeling was executed.

## Disabled entry and integration corrections

Original inter-image ionic gap is16.00026725 A. Expanded cell adds10 A along z and translates every ion+5 A, yielding26.00026725 A gap; in-plane cell, element order, IDs, markings and relative ionic positions remain fixed. Lower/upper margins and full records are in input_manifest.json. Pairwise relative geometry tolerance1e-12 A; extxyz roundtrip tolerance1e-8 A. No relaxation.

Six label-free proposals comprise original/expanded vacuum crossed with: all-periodic baseline; slab pbc=(True,True,False) without correction; same slab with correction. Use periodic baseline vs periodic expansion for vacuum sensitivity; same-cell slab-off vs slab-on for dipole sensitivity; periodic vs slab-off diagnoses boundary effects separately. Existing 0-vacuum proposals are copies for a review pack, not authorization to repeat completed baselines.

PW-compatible API for GPAW26.7: `poissonsolver={'dipolelayer':'xy'}` with pbc_z=False and z perpendicular to the in-plane vectors. Evidence: `old/pw/hamiltonian.py` constructs DipoleCorrection and calls check_direction; `dipole_correction.py` rejects periodic correction direction. New `new/pw/builder.py` forwards solver parameters and `new/pw/poisson.py:DipoleLayerPWPoissonSolver` requires exactly two periodic directions and chooses the remaining axis. The truthy layer argument must correspond to xy in this z-normal geometry. Do not reuse the historical all-PBC verifier unchanged.

All records have launch_enabled=false, no owner, unresolved mesh/width, PW500/PBE preserved. No executable launch runner is included. A must integrate explicit per-record width/mesh/PBC/Poisson metadata, full source hashes and output verifiers, exact owner/label registration, immutable/no-existing-run guards, four-rank resource and disk guards; verify exact relative geometry and cell change before any launch. Actual DFT remains blocked pending the fixed-width numerical decision after A6/8. Only static preparation and read-only numerical audits ran.

Reproduce from /workspace/-:

```bash
/workspace/.venvs/gpaw-mpi/bin/python research/mace-v12-transfer/coordination/reports/window-b/v19_numerical_target_and_vacuum_entry_review/audit.py
```
