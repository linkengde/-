# Stage C-D6 — Minimal 77-atom Ag/Ti₃SiC₂ DFT input preparation

## Status

```text
DFT_INPUT_READY=YES
DFT_EXECUTION_READY=NO
DFT_STARTED=NO
LAMMPS_STARTED=NO
TSCBOP_OR_MORSE_MODIFIED=NO
FITTING_STARTED=NO
```

The three fixed-geometry candidate structures were copied from the D5 geometry-audited archive and checked against its SHA-256 manifest. The copy is byte-identical; no coordinates, atom order, cell, PBC, composition, or interface distance were changed. Each candidate has a JSON input specification and is bound to `run_fixed_singlepoint.py`. All configurations explicitly contain `execution_authorized=false` and `launch_enabled=false`; there is no execution authorization in this package. The runner refuses to construct GPAW unless a separate signed release binds the actual config, structure, and runner hashes.

`DFT_INPUT_READY=YES` means the structural and method input files are assembled for review. It does not mean the calculations are released, that the chosen k mesh is converged, or that the present machine has sufficient resources. The current resource record remains `UNVERIFIED`; no DFT, LAMMPS, MD, TTM, fitting, or potential edit was performed.

## Structures and fixed-geometry scope

All three structures have 77 atoms, composition Ag₂₈Ti₂₁Si₁₄C₁₄, a common cell, PBC `T T T`, and atom IDs 1–77. The substrate is the seven-plane Si–Ti–C–Ti–C–Ti–Si slab; the Ag film is four ABC-A-stacked Ag(111) planes. The coherent interface is the D5 0° shared-strain model; it is **not** the archived 12.51983° low-mismatch interface. The complete lattice vectors, per-structure file/config hashes, gaps, and strain values are in [`STRUCTURE_LIST.csv`](STRUCTURE_LIST.csv). The input geometries are unrelaxed.

| Task ID | Structure | Fixed plane gap | Min Ag–Si | Min Ag–Ti | Min Ag–C | Top vacuum |
|---|---|---:|---:|---:|---:|---:|
| `C-D6-M77-GAP-2P16A` | `Ag_Ti3SiC2_M77_gap_2p16A.extxyz` | 2.16 Å | 2.160 Å | 4.530 Å | 5.585 Å | 19.840 Å |
| `C-D6-M77-GAP-2P40A` | `Ag_Ti3SiC2_M77_gap_2p40A.extxyz` | 2.40 Å | 2.400 Å | 4.753 Å | 5.814 Å | 19.600 Å |
| `C-D6-M77-GAP-2P64A` | `Ag_Ti3SiC2_M77_gap_2p64A.extxyz` | 2.64 Å | 2.640 Å | 4.977 Å | 6.043 Å | 19.360 Å |

The common in-plane cell vectors are `(7.4468326699, −2.5796585079, 0) Å` and `(−1.4893665340, 7.7389755238, 0) Å`; the third vector is `(0, 0, 37.9451596) Å`. The cell volume is 2041.024552 Å³. The applied biaxial strains relative to the frozen D5 references are Ti₃SiC₂ −3.00446% and Ag +3.09752%. The 2.40 Å separation is a geometric reference, not a relaxed equilibrium distance. The three geometries differ only by rigid z translation of the complete Ag film, preserving its internal coordinates and layer spacing.

The D5 geometry-only audit found a minimum all-atom separation of 2.053774 Å and no pair below 1.5 Å. Its separate LAMMPS force admission gate was not measured and remained HOLD; this DFT input package does not convert that gate into a pass.

## Candidate GPAW input method

The three JSON configs use one common candidate method so their fixed-geometry energy/force trends are comparable:

- GPAW 26.7.0, ASE 3.29.0, `gpaw-data` 1.2.1; PBE with PAW.
- Plane-wave cutoff 500 eV; Fermi–Dirac width 0.10 eV.
- SCF energy, density, and eigenstate tolerances each `1e-5`; maximum 160 iterations.
- Pulay mixer: beta 0.05, `nmaxold=8`, weight 100.
- Proposed first mesh: 4×4×1, the D5 resource plan's comparable-density starting proposal; GPAW default integer Monkhorst–Pack convention, no explicit Gamma override or extra offset.
- Four MPI ranks and one OpenMP thread.
- Fixed ions and fixed cell; no relaxation or cell change.
- Requested outputs: native total energy, force-consistent free energy, complete atomic force vectors, SCF status/iterations, and full stress tensor. Slab stress is vacuum-diluted and is diagnostic only unless a separate convention is approved.

The 4×4×1 mesh has **not** been convergence-tested for these M77 structures. The three-gap series must use the same mesh and method. A denser central-gap comparison (5×5×1) is a possible later convergence check, not part of an authorized calculation in this task. Do not interpret the three gap points as TRAIN, DEV, or TEST labels; their data role remains unassigned until a separate review.

## Intended scientific use and limitations

If later calculated, converged, and audited, the fixed-geometry series can inform Ag/Ti₃SiC₂ interfacial interaction trends and a cross-potential trend check across the three prescribed separations. It is a compact, strained, coherent 0° model. It cannot by itself validate the 12.51983° low-mismatch interface, rank independently relaxed registries, establish a production cross potential, or prove transferability to thermal/liquid interfaces. Energies from different model families or different numerical settings must not be compared as if their energy zeros and methods were identical.

## Package integrity and release boundary

`DFT_INPUTS/structures/` contains the unchanged D5 extxyz files; `DFT_INPUTS/configs/` contains one explicit non-authorizing method config per structure. `INPUT_MANIFEST.json` records the source manifest, structure/config hashes, method, and runner hash. `RESOURCE_ESTIMATE.md` gives a relative planning proxy only. `SHA256SUMS.txt` covers every package file except itself and this checksum file.

Before any future execution, a separate review must approve the k-point convention/convergence plan, fresh resource and storage preflight, exact setup hashes, and signed execution authorization bound to the final config/structure/runner hashes. The current machine's 32-GiB history and unknown peak-memory requirement are not cleared by preparing these inputs.
