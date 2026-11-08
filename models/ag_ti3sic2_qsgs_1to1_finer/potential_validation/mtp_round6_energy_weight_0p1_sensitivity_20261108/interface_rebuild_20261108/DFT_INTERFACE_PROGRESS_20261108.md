# Ag(111)/Ti₃SiC₂(0001) DFT interface progress

Date: 2026-11-08. This note distinguishes prepared starting geometries and preliminary SCF checks from a validated interface separation curve. No binding energy, work of separation, or DFT-relaxed interface value has yet been accepted.

## Interface cells prepared

The geometric pre-screen used the repository's existing PBE reference cells (`a_Ag = 4.152 Å`; `a_MAX = 3.068 Å`, `c_MAX = 17.669 Å`). These values have not yet been re-optimized with the final GPAW LCAO settings, so the match is provisional. The best small coincidence candidate is 13 Ag(111) surface sites matched to 12 Ti₃SiC₂(0001) basal cells, common surface cell about 10.6067 Å, area 97.43 Å², with balanced strain about 0.200%. The 6-plane MAX slab plus 3-layer Ag film contains 111 atoms. A 12-plane thickness check contains 183 atoms. The 1×1 method-control cell has about 2.25% balanced strain and is not the physical low-strain candidate.

Prepared inputs cover top C/Ti/Si termination, 3/13, 6/13 and 13/13 Ag-layer coverage, and top/two hollow registries. Eleven fixed-composition separation paths use nominal gaps 2.2, 2.6, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0 and 10.0 Å. All are unrelaxed starts. The selected preliminary primary path is the low-strain 13:12 Ti-terminated, top-registry 6-plane slab.

This record preserves the cell dimensions, strain, coverage/registry matrix, and scan plan. The geometry generator and candidate manifest are local working files and are not included in this documentation-only branch.

## Preliminary DFT checks

The GPAW method is PBE/LCAO/dzp, non-spin-polarized, Fermi–Dirac width 0.1 eV, Pulay mixer (`beta=0.03`, `nmaxold=8`, `weight=100`), with energy/density/eigenstate thresholds 10⁻⁵. The current sparse actual-interface diagnostic uses Γ only; it is not a converged k-point-quality production value.

| Check | Result | Interpretation |
| --- | --- | --- |
| 1×1 Ti/top method-control path, 3×3×1 | Nominal gaps 10, 8, 6 and 5 Å each converged; the 5 Å point took 44 SCF cycles | Tests the code and electronic settings only. This cell has 2.25% strain and only 9 atoms; it is not the low-mismatch interface. |
| 13:12 low-strain Ti/top, 3×3×1 at 10 Å | Stopped after 10 minutes before the first SCF iteration; no energy/force accepted | 2,565 DZP basis functions, 748 bands and 5 irreducible k-points made this initial screen too costly in the local run. It is not evidence of non-convergence. |
| 13:12 low-strain Ti/top, Γ only at 10 Å | SCF converged in 44 cycles: native energy −576.23719063 eV/cell; free energy −578.34881261 eV/cell; Ag-film total Fz −0.02559 eV/Å; maximum atomic force 1.2884 eV/Å | A valid Γ-point SCF diagnostic on the unrelaxed starting geometry. The large residual force means this is not a relaxed surface or interface reference. Do not use it to calculate binding energy. |
| 13:12 low-strain Ti/top, Γ only at 3 Å | Resumed after the saved 10 Å point. At the time of this repository update (05:13 UTC), 29 SCF iterations had completed; the latest electronic energy was −587.732560 eV/cell, but the point was not converged. | In progress; no 3 Å energy or force is accepted until GPAW reports convergence. |

The 10 Å point was written before the explicit post-save MPICH barrier failed with `Message truncated`. The continuation therefore skips that per-point barrier. The successful launch uses the MPICH `mpiexec.hydra` shipped with the GPAW environment; the system OpenMPI launcher is ABI-incompatible with this GPAW build. The failed 3×3 and barrier attempts are logged separately; no failed or unfinished energy was copied into the accepted point list.

The run logs and progress JSON files remain in the local working directory and are not included in this documentation-only branch. The 3 Å Γ-point continuation remained active at 05:13 UTC; until GPAW reports SCF convergence, the partial energy is not an accepted contact energy and must not enter a binding-energy calculation.

## What remains before an interface curve is admissible

The 13:12 candidate needs bulk lattice checks under the final interface settings, surface-thickness and k-point convergence, ionic relaxation at each constrained separation, and matched separated Ag-slab/MAX-slab reference energies under identical strain, atom counts, cell and settings. The actual post-relaxation gap and forces must be recorded. If a periodic construction contains more than one interface, the interface count must be included when converting energy to area-normalized separation work.

The local scanner is deliberately a rigid static SCF screen; it does not optimize ions. No DFT interface binding curve or W_sep has been obtained. The generator, scanner, and their run files are not part of this summary-only upload. No energy-and-force MTP refit, new heating MD, or two-contact calculation is authorized by this stage of the workflow; those remain behind the DFT interface gates.
