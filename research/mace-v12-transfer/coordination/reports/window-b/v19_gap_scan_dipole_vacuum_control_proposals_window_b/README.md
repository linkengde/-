# AgSi dipole and vacuum control geometry proposals

These are label-free, launch-disabled input proposals derived from the registered AgSi parent gap-scan geometry. They do not authorize or start a calculation. The registered gap-scan manifest and its source input are unchanged.

`AgSi_26A_vacuum_parent_dipole_OFF_PROPOSAL.extxyz` is a byte-for-byte copy of the registered parent geometry. Its proposed method changes only `poissonsolver` from `{"dipolelayer":"xy"}` to `{}`.

`AgSi_36A_vacuum_parent_dipole_ON_PROPOSAL.extxyz` preserves all atom identities and relative coordinates, increases only the z cell vector by 10 Å, and shifts every atom +5 Å in z. It keeps the registered dipole-ON method. This changes the geometric empty length `cell_z − (z_max − z_min)` from 26.000267 Å to 36.000267 Å while keeping the slab centered. This coordinate-span measure is not an electron-density-defined vacuum width.

Both proposals retain 52 atoms, composition, atom order, IDs, pair markers, in-plane cell, PBC `[true,true,false]`, and the shared uploaded-CIF parent hash. They are numerical boundary controls from the same parent, not independent validation structures. The extxyz files contain no energies or forces. Both proposal records have `launch_enabled=false`, `owner=null`, and a diagnostic-only role.

`prepare_control_proposals.py` verifies the registered source hashes and creates the two files without overwriting existing outputs. `proposal_manifest.json` records exact method and source/proposal hashes; `verification.json` records geometry and label-free checks. To reproduce from the repository root:

```bash
python3 research/mace-v12-transfer/coordination/reports/window-b/v19_gap_scan_dipole_vacuum_control_proposals_window_b/prepare_control_proposals.py
```

No DFT calculator, MPI, or GPAW import is used.
