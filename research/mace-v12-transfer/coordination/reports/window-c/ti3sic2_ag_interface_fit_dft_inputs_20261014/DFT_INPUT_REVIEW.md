# DFT input and release review

## Geometry inventory

`INPUT_MANIFEST.json` and `STRUCTURE_FAMILY_LABELS.csv` list 15 proposed fixed geometries: 6 bulk Ti₃SiC₂ and 9 Ag/Ti₃SiC₂ interface gaps. `build_input_package.py` checks the E7 and 217-atom source hashes, stoichiometries, atom counts, and extxyz round-trip identity. It does not evaluate energy or forces.

For the interface scan, the slab cell, PBC, lateral registry, atom order, atom IDs and in-plane strain are retained. Only all Ag z coordinates are translated by `target_gap - 2.0 Å`. No atom or cell relaxation, lateral registry shift, or geometry optimization occurs. Cross-interface distance metrics are recorded in each config.

## DFT input contents

Each `configs/*.json` binds the exact geometry SHA, family, task ID, perturbation, method and requested labels. Both systems request four MPI ranks and one OpenMP thread per rank. The configs bind the audited PBE PAW setup hashes: Ti `bf822f…add2f`, Si `2defcf…54ef0`, C `fbef4c…3c983`, and Ag `bff00a…145e97` (Ag is needed only for the interface). `run_fixed_singlepoint.py` checks ASE/GPAW/gpaw-data versions, `OMP_NUM_THREADS=1`, PAW setup hashes, and an external authorization binding task ID plus config, structure and runner hashes before calculator construction. The package has no authorization file and every config sets `execution_authorized=false`.

Requested labels:

- Bulk Ti₃SiC₂: native total energy, force-consistent free energy, complete per-atom force vectors, and the complete ASE stress tensor in eV/Å³.
- Ag/Ti₃SiC₂ slab: the same energy and force outputs; full stress is recorded only as a vacuum-diluted diagnostic and must not be mixed with bulk stress in a fit loss.

Native energy (`force_consistent=False`) and free energy are distinct fields and must not be interchanged. A non-converged calculation must not be ingested as a training label.

## Resource and method gates

- The Ti₃SiC₂ settings match the E6/B archived bulk method and mesh. The 100-atom slab settings in E6 had not completed their own convergence review; this package does not use those slab inputs.
- The 217-atom interface's 2×2×1 Gamma-centered mesh is inherited from a prior proposal, not from a demonstrated k-point convergence test. Before a full scan, compare at least one identical fixed geometry against a denser in-plane mesh and agree an error budget.
- The previous 394-atom DFT attempt exceeded a 32 GiB memory limit during initialization. The 217-atom model was selected as a smaller candidate, but no matching successful memory benchmark exists. Resource availability and wall time are therefore `UNVERIFIED`; do not launch the interface batch on the basis of the older 394-atom failure or the smaller atom count alone.
- No exact cost estimate is assigned to the 217-atom scan without a successful benchmark. Existing 48-atom runs imply roughly 3–4.2 hours per bulk point on four ranks; six points are a planning range of about 18–30 rank-hours, not a guarantee. The nine interface points are not budget-released.

## Release state

```text
STRUCTURE_INPUTS=PREPARED
DFT_LABELS=NOT_GENERATED
DFT_EXECUTION_AUTHORIZED=NO
INTERFACE_KPOINT_CONVERGENCE=NOT_ESTABLISHED
INTERFACE_RESOURCE_BUDGET=UNVERIFIED
TSCBOP_MODIFIED=NO
FITTING_STARTED=NO
DATASET_FIT_READY=NO (no labels and no independent DEV/TEST family)
```
