# DFT input preparation: Ti₃SiC₂ backbone and Ag/Ti₃SiC₂ interface

## Result

Prepared **15 fixed-geometry DFT input candidates**: six Ti₃SiC₂ bulk response points and nine Ag/Ti₃SiC₂ Si-terminated interface-distance points. Geometry validation passes for all 15. This is an input package, not a completed DFT dataset: **no new DFT energy, force, or stress labels exist yet**. Every config has `execution_authorized=false`; no GPAW calculation, TSCBOP edit, or fit was performed.

The six bulk geometries and nine interface geometries are under `inputs/`; the per-task parameter locks are under `configs/`. The family and role ledger is `STRUCTURE_FAMILY_LABELS.csv`; the split policy is `TRAIN_DEV_TEST_SPLIT_PLAN.md`; provenance and SHA bindings are in `INPUT_MANIFEST.json` and `SHA256SUMS.txt`.

## Prepared coverage

| Priority | Candidate set | Count | Parent family | Intended labels | Status |
|---|---|---:|---|---|---|
| 1 | Volume −2%, +2%; shear −0.5%, +0.5% | 4 | `TSC_KISI_COD9009647_BASELINE_DERIVED_FAMILY` | Native energy, free energy, forces, full stress | Geometry only; TRAIN candidate pending role freeze |
| 2 | Si atom `local_new_id=7`, z ±0.02 Å | 2 | Same TSC baseline-derived family | Native energy, free energy, forces, full stress | Geometry only; TRAIN candidate pending role freeze |
| 3 | Ag slab rigid z scan, Si-terminated interface plane gap 1.8–3.4 Å at 0.2 Å increments | 9 | `AG_TSC_SI_TERM_TWIST12P51983_217_PARENT_FAMILY` | Native energy, free energy, all forces; slab stress diagnostic only | Geometry only; future unified-model TRAIN candidate, method/resource hold |

The interface scan uses the archived 217-atom structure unchanged as its parent: Ag84Ti57Si38C38, 12.51983° twist, +1.1029% Ag in-plane strain, 0% Ti₃SiC₂ strain, 20 Å upper vacuum, TTT PBC. The 1.8–3.4 Å values are **plane gaps**, not nearest-pair distances or optimized `d0`. For each child, only the z coordinates of all Ag atoms change by `gap−2.0 Å`; cell, PBC, lateral registry, atom order/IDs, and every Ti/Si/C coordinate remain identical. The minimum Ag–Si distance spans 1.8–3.4 Å; all checked Ag–Ti/C contacts remain larger, and no periodic pair is below 1.5 Å.

## DFT settings

The Ti₃SiC₂ inputs match the audited E6/B bulk setup: GPAW 26.7.0, ASE 3.29.0, gpaw-data 1.2.1, PBE/PAW, PW 500 eV, Fermi–Dirac 0.10 eV, Pulay β=0.05 (`nmaxold=8`, `weight=100`), energy/density/eigenstate thresholds 1e-5, maxiter 160, and the 10×10×4 Monkhorst–Pack mesh. Fixed geometry; native energy and force-consistent free energy remain separate; full forces and stress are requested.

The interface candidates use the earlier proposed 2×2×1 Gamma-centered mesh with the same functional, cutoff, occupations, mixer and thresholds. This mesh has no convergence evidence for this cell. It must be tested against a denser in-plane mesh and released with an agreed error budget before any production scan. Full 3D slab stress is reported only as a vacuum-diluted diagnostic, not mixed into the bulk stress loss.

## Family-level split

- **TRAIN candidates:** keep the complete baseline-derived bulk family together. The six new geometries are candidates only; current E5 baseline/Ti/C force labels remain numerical controls and are not silently relabeled or copied. Any role migration requires a separate freeze.
- **Future unified-model TRAIN candidate:** keep all nine interface gaps together as one family. These Ag-containing structures cannot constrain the three-element Ti–Si–C TSCBOP. They are for a separately scoped unified Ag–Ti–Si–C model.
- **DEV:** no qualifying independent bulk or complete interface source family is currently available. Leave `SOURCE_PENDING`; do not split children of either proposed TRAIN parent into DEV.
- **TEST:** leave evaluator-controlled. No sealed structure or label was read or enumerated.

Therefore there is no populated, independent DEV or TEST set. These inputs cannot yet support a controlled parameter-selection experiment or a generalization claim.

## Resource and release status

The prior 48-atom matched bulk calculations took about 2.83–4.20 hours per point on four MPI ranks. A planning estimate for six points is 18–30 rank-hours, not a guarantee. No successful resource benchmark exists for this 217-atom interface model. A prior 394-atom interface single point exceeded a 32 GiB memory environment; that does not prove the 217-atom case fits. Interface runtime, memory and checkpoint needs remain unverified.

`run_fixed_singlepoint.py` requires an external authorization JSON binding task ID, config SHA, structure SHA and runner SHA. No authorization is included. Its no-authorization refusal was tested; the script has not been used to start DFT. The interface method/resource review, independent DEV/TEST sources, explicit TRAIN role freeze and actual DFT execution authorization remain open.

## Final state

```text
DFT_INPUTS_PREPARED=YES
GEOMETRY_VALIDATION=PASS (15/15)
DFT_LABELS_GENERATED=NO
TRAIN_ROLE_FROZEN=NO
INDEPENDENT_DEV=NO
SEALED_TEST_ACCESSED=NO
INTERFACE_KPOINT_CONVERGENCE=HOLD
DFT_EXECUTION_AUTHORIZED=NO
TSCBOP_MODIFIED=NO
FITTING_STARTED=NO
```
