# Stage67 stable meltpool liquid Ag extraction report

## Status

```text
MAXIMUM_CLUSTER_EXPORT = PASS
FULL_LIQUID_AG_EXTRACTION = BLOCKED
OVERALL_EXTRACTION_STATUS = PARTIAL
LAMMPS_RUN = NO
SOURCE_DUMP_MODIFIED = NO
```

## Source and frame

- Source: Stage67, 200 Å transverse model, audited stable-meltpool reference branch.
- Audited window: 0.155–0.165 ps (step 1550–1650; 11 sampled states at 1 fs intervals). This is the existing 10 fs window assessment, not a claim of long-time or thermodynamic stability.
- Extracted raw frame: step 1650 / 0.165 ps.
- Cell: 150.1 × 200 × 200 Å; 303,372 atoms.
- Boundary flags: `pp pp pp` (periodic in all three directions).
- Source frame SHA256: `91D04A6FA3CB1B2FD9EB668B392CAD04ECB593620723BEAC3E12EA4E7B8E85FC`.

## Counts and cluster result

- Total Ag atoms (`type=4`), counted from the source dump: **94,624**.
- Existing step1650 Q6/PTM summary, `N_LiquidLikeAg`: **8,728**.
- Existing maximum liquid-like Ag cluster: **1,032 atoms**.
- Maximum-cluster center of mass from the extracted coordinates: **(68.920100, 88.326150, 89.191421) Å**.
- Existing maximum-cluster bounds: x = 50.477816–99.775597 Å; y = 60.367842–119.957271 Å; z = 60.215459–139.711422 Å.
- Maximum-cluster `PBC_contact`: **NO**.
- Maximum-cluster `PBC_spanning`: **NO**.

These PBC flags describe the audited maximum cluster only. The simulation cell itself is periodic in x, y, and z.

## Frozen classification provenance

The existing Stage67 implementation selects Ag (`type=4`), applies OVITO Polyhedral Template Matching to the selected Ag atoms, and forms the candidate mask from `Structure Type == 0` (PTM `Other`), `CellTl_valid == 1`, and `CellTl >= 975 K`. The existing summary records `Q6_threshold=0.446766` and `cluster_cutoff_A=3.74625`; in the available implementation, the numeric Q6 threshold is summary metadata while the per-atom mask uses the PTM `Other` proxy. No criterion was changed here.

The saved source package contains the aggregate count and IDs for the maximum 1,032-atom cluster, but does not contain per-atom PTM/Q6 labels or IDs for the remaining liquid-like Ag atoms. The dump itself contains coordinates, velocities, `CellTl`, and `CellTl_valid`, but no Q6/PTM column. The available OVITO Python bindings are built for CPython 3.14; no compatible interpreter or `ovitos.exe` was found. Therefore the full 8,728-atom selection could not be recomputed and verified against the existing count.

## Exported files

- `meltpool_cluster_step1650.dump`: 1,032 Ag atoms from the existing maximum-cluster ID list; preserves source atom records and box bounds; directly loadable as a LAMMPS dump in OVITO.
- `meltpool_cluster_step1650.xyz`: the same 1,032 Ag atoms in extended XYZ with atom IDs, step/time, and periodic lattice metadata; loadable in OVITO.
- `liquid_Ag_step1650.dump`: **not generated**. It must represent all 8,728 classified liquid-like Ag atoms; substituting only the 1,032-atom maximum cluster would be incomplete and misleading.

## Data integrity and actions

- Source dump opened read-only; its SHA256 matches the existing input-package manifest.
- The 1,032 IDs were matched against the source frame; all matched atoms are type 4 Ag.
- Extracted coordinate center agrees with the existing summary cluster COM to rounding precision.
- No LAMMPS/MPI was started; no original dump, restart, potential, or input was modified.
