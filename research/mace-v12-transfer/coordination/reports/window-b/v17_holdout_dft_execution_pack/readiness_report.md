# Readiness report: V17 holdout DFT execution pack

## Status

Preparation is complete; **no DFT process, GPAW import, MPI launch, MACE inference/training, MD or TTM was run**. No authorization file was created. All nine prospective labels remain blocked pending A's role review and explicit owner assignment.

Disclosure: while locating the pinned source settings, an early broad recursive text search also matched some existing V17 calculation log/summary text and emitted result lines in command output. I stopped the search immediately, did not copy or use those result values in this pack, and did not read any further DFT result files. This was an accidental deviation from the assignment's result-read boundary; this report preserves it for A's review.

## Geometry and provenance checks

- Verified the source `final_set/SHA256SUMS.txt` inventory; all listed candidate bundle files passed.
- Copied the nine candidate extxyz files without rewriting. The package builder checks that each copied byte stream exactly matches its source-file SHA256.
- `prepare_pack.py` checks atom count, formula, ordered elements and LAMMPS IDs, cell, PBC, and absence of energy/force properties. It also hashes the geometry-only direct and parent input files against the candidate provenance records.
- The manifest retains the earlier geometry-audit decision and species-specific contact/duplicate evidence. Historical labeled split files were not reopened; the accidental search described above did match existing calculation-output text.
- The upstream geometry audit reports nine geometries with no exact/near inventory hits and no role-set overlap. Its audit reports the compared inventory as 197 frames across 18 model-split files and 38 published inputs. Candidate motifs are inherited from existing small-cluster parents; morphology and thermal-state independence are not established.

## Proposed ownership and label policy

Proposed only: A owns the three development-validation calculations; B owns the six withheld-test calculations. A must explicitly accept or revise each assignment. For withheld-test labels, A must publish a hashed frozen-checkpoint/selection record before reading the labels. B must not analyze, compare, or model-score withheld-test values; its future archive role is mechanical identity/convergence/finite-value verification only.

## Guards prepared

The future runner requires explicit paths for workspace, pack, inputs, scratch run root, archive root, authorization, owner and label list. It checks the A-issued authorization against the prepared role-manifest and exact input hashes; enforces role/owner assignment; checks geometry identity, IDs/order, cell/PBC and absence of labels; refuses an existing per-label output or archive and stale staging paths; takes a nonblocking execution lock; refuses concurrent MPI/DFT; checks four actual CPU cores and GPAW 26.7.0 MPI/ScaLAPACK features; sets one OpenMP/BLAS thread per MPI rank; and runs labels serially with four ranks.

The GPAW draft retains the V17 PBE/PW500/Gamma/FermiDirac0.1 recipe, existing Pulay mixer and convergence thresholds, fixed coordinates, native extrapolated energy plus separate free energy when available, progress/log evidence and a local-only `state.gpw`. The archive draft requires convergence evidence, finite energy/forces, exact source hash and structural identity, and archives only compact extxyz, summary, progress, GPAW/launcher logs, verification and hashes.

## Static checks

The system `python3` lacked ASE, so geometry-only manifest preparation used the existing `/workspace/.venvs/gpaw-mpi/bin/python` (Python 3.12.14, ASE 3.29.0). The builder imported ASE only; GPAW was not imported. The candidate geometry manifest was built from the nine source inputs and their geometry-only provenance. Shell parsing and Python AST parsing are used only to check draft syntax; the guarded runner, GPAW driver and archive verifier are not executed. No labels are present in the prepared inputs or manifest.

## Remaining gate

Do not launch until A reviews/finalizes the roles and issues a separate exact-hash authorization/owner assignment. This pack is not a training authorization or a claim that these candidates are independent of their small-cluster parents.
