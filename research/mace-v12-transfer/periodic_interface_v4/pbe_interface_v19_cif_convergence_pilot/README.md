# CIF-supported matched-kmesh numerical pilot

A approved the same52atom AgSi fixed-ion registry for Gamma(A) versus2x2x1(B) comparison after B source and entry reviews. This is a numerical probe, excluded from training/validation integration. Both jobs require exact registered owner identity. Source CIF is attributed COD9009647/Kisi et al1998; remote checksum not independently retrieved.

Queue checks launch/owner/input before claiming; direct calculator also enforces owner registration. Existing archives are immutable and their complete inventories verified; archive plus local run coexistence stops for same-job review. Full method/version/energy convention/progress/force-summary checks enforced. Local state.gpw remains in workspace run directory, excluded from compact publication. No automatic checkpoint deletion. Startup requires3GiB free; SCF hook checks1.5GiB reserve before checkpoint writes.

Compare native/free energy, all-atom force differences and marked-pair projection. Two meshes do not establish full convergence. Vacuum/dipole/thickness/Ag-strain effects and common-target consistency remain follow-ups. Fixed-ion stress probes are not equilibrium surface evidence or MACE accuracy evidence.

A additional3x3x1 pilot is registered for the same frozen geometry. When more than one owner record exists, specify exact label to run_queue.py; the completed Gamma job is not relaunched. compare_mesh_results.py accepts two exact label arguments and writes a separate pairwise file. CompareGamma/2x2 and2x2/3x3; vacuum/dipole convergence remains separate.
