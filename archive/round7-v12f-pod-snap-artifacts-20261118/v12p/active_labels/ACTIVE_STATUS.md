# v12o hotspot labels and v12p continuation

Status snapshot: 2026-11-19 07:44 UTC. All four GPAW PBE/PW350 Gamma single-point labels for local crops from the 458,248-atom SNAP v12o no-laser 10 fs endpoint are SCF-converged. Settings were density/eigenstates/energy 1e-5, maxiter 200.

| Crop | SCF iterations | Atoms | DFT force targets |
|---|---:|---:|---:|
| Ag-Si 430045/250182 | 52 | 44 | 2 |
| Ag-Si 455657/339883 | 40 | 38 | 2 |
| Ag-C 454089/336261 | 37 | 40 | 2 |
| Ag-Ti 412351/186486 | 42 | 37 | 2 |

All four JSON records say `CONVERGED`; the logs contain `Converged in ...`; extxyz target IDs and finite force arrays passed preparation checks. The candidate preserves the 128 train / 12 validation split and keeps Si01 center atom 104359 held out. It has 30 local force frames, 158 training configurations total, and 33 targeted DFT force rows.

Pipeline session 52808 is now fitting v12p POD/SNAP. Log: `../v12_candidate_forceorder_fixed/v12p_pipeline.log`. After fitting it evaluates Si01, runs full-model static force checks, and gates no-laser 10 fs cold screens at 500 eV/Å / finite forces. Heating remains blocked pending stability.
