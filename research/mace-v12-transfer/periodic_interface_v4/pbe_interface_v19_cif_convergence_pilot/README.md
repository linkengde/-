# CIF-supported matched-kmesh numerical pilot

A approved the same52atom AgSi fixed-ion registry for Gamma(A) versus2x2x1(B) comparison after B source and entry reviews. This is a numerical probe, excluded from training/validation integration. Both jobs require exact registered owner identity. Source CIF is attributed COD9009647/Kisi et al1998; remote checksum not independently retrieved.

Queue checks launch/owner/input before claiming; direct calculator also enforces owner registration. Existing archives are immutable and their complete inventories verified; archive plus local run coexistence stops for same-job review. Full method/version/energy convention/progress/force-summary checks enforced. Local state.gpw remains in workspace run directory, excluded from compact publication. No automatic checkpoint deletion. Startup requires3GiB free; SCF hook checks1.5GiB reserve before checkpoint writes.

Compare native/free energy, all-atom force differences and marked-pair projection. Two meshes do not establish full convergence. Vacuum/dipole/thickness/Ag-strain effects and common-target consistency remain follow-ups. Fixed-ion stress probes are not equilibrium surface evidence or MACE accuracy evidence.

A additional3x3x1 pilot is registered for the same frozen geometry. When more than one owner record exists, specify exact label to run_queue.py; the completed Gamma job is not relaunched. compare_mesh_results.py accepts two exact label arguments and writes a separate pairwise file. CompareGamma/2x2 and2x2/3x3; vacuum/dipole convergence remains separate.

Gamma-to2x2 pilot exceeds all provisional budgets. B4x4 follow-up is registered; compare2x2/4x4 and3x3/4x4 to distinguish grid-density and even/odd sampling sensitivity. These comparisons establish only this one geometry/settings case; no transfer of convergence claim to other surfaces or compact training frames.

## New independent owners: 5x5/6x6 continuation

Both fresh owners are registered. See `new_owner_k5_k6_registration.json` for pinned owners, source, method and unchanged execution-file hashes. Old four records and their owners are preserved. A runs only `AgSi_COD9009647_pilot_k5x5`; B runs only `AgSi_COD9009647_pilot_k6x6`, each in its own environment:

```bash
cd /workspace/-
. /workspace/.setup/activate-gpaw.sh
python3 research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_cif_convergence_pilot/run_queue.py window-a AgSi_COD9009647_pilot_k5x5
# B uses window-b AgSi_COD9009647_pilot_k6x6 instead.
```

After verified publication compare 4/5, 4/6 and 5/6 with `compare_mesh_results.py` using exact labels. Keep existing numerical budgets. These labels remain numerical-only; no V19 training. Check smearing if mesh differences still exceed budgets, then vacuum/dipole separately.
