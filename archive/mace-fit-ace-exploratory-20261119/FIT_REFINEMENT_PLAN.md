# MACE-derived ACE/MTP fit: next iteration

This is a residual-driven plan for the MACE-fit data candidate. It does not use Round7 labels and does not replace the baseline ACE or MTP models.

## What the residuals show

- ACE passes the aggregate development RMSE gates: 5.55 meV/atom energy and 0.0930 eV/Å force-component RMSE. It does not pass every frame individually: the v25 `rms0p045` and `rms0p06` bulk displacement frames score 0.156 and 0.199 eV/Å force-component RMSE.
- In the `rms0p045` frame, the per-species component RMSE is 0.194 eV/Å for Ti and 0.125 eV/Å for C. In `rms0p06`, C is 0.274 eV/Å and Ti is 0.167 eV/Å. Si remains lower in those two frames. This points to missing force-response coverage for displaced Ti–C environments, rather than a general Si-force problem.
- The `Ti_gap_2p3_rank_recovery_v17_01` training frame itself has reference-force vector RMS 3.10 eV/Å (component RMS 1.79 eV/Å) and maximum atomic force 5.68 eV/Å; its baseline ACE force-component RMSE is 0.831 eV/Å. Its Gamma SCF converged and archive hashes pass, but the source geometry declared a 1×4×2 mesh, so k-point convergence remains unresolved. Do not discard it solely because its residual is high.
- The exact Si hotspot-center gate is still unscored because the MACE-fit split does not register hotspot-center atom IDs. The maximum error across all Si sites (0.245 eV/Å on development) is only a proxy.

## Recommended next actions

1. Freeze the original ACE/MTP files and the seven-frame development set as the baseline.
2. The provenance audit is complete: the V17 `Ti_gap_2p3` Gamma label is SCF-converged and hash-verified, but its declared 1×4×2 source mesh leaves k-point convergence unresolved. A 1×4×2 rerun was stopped after two SCF iterations because it was too costly in this runtime; no result from it is used.
3. A sensitivity fit that scales only this frame's force weight to 0.05 meets the stated aggregate train/development gates and every current development frame's energy/force limits. Its own `Ti_gap` training force RMSE remains 0.914 eV/Å, so this is suitable only if this rare high-force geometry is outside the target dynamics. The `.ace.yaml` and `.yace` passed a small LAMMPS compatibility check. Keep this as an alternate exploratory candidate, not a production-qualified model.
4. Do not add the proposed 12–20 displaced-bulk DFT structures yet. First decide whether near-gap, high-force `Ti_gap` environments are part of the intended simulation domain. If they are, obtain a consistent, verified reference for this structure and several neighboring configurations; if they are not, retain the downweighted model for exploration and focus new labels on the actual target temperature/impact regime.
5. Do not transfer the ACE 0.05 force downweight to MTP: an MTP10 fit that removed the row scored 78.5 meV/atom and 0.495 eV/Å on development. Keep the full-data MTP10 as the baseline while diagnosing architecture/coverage; its 0.269 eV/Å force RMSE still fails.
6. If this development set continues guiding model changes, reserve a new provenance-grouped independent test set before claiming final qualification.
7. Register the intended Si hotspot-center atom IDs and score the requested center-vector metric explicitly.
8. The million-atom LAMMPS static check is complete (1,021,793 atoms, 4 MPI ranks, 35 s for `run 0`). The force audit finds 330,452 atoms above 0.5 eV/Å, while the minimum nearest-neighbor distance is 2.088 Å and no pair is below 1.5 Å. Before dynamics, extract representative high-force Ti–C and Ag–Ti/C/Si local environments, obtain reference labels, and check local minimization stability. Then validate smaller cells for solid/liquid Ag, interfaces, and elevated-temperature behavior; measure million-atom MD throughput only after those checks pass.

The current development energy RMSE passes the aggregate target, but one v19 frame has an absolute energy error of 10.86 meV/atom. Thus per-frame energy compliance also needs to be distinguished from the stated aggregate RMSE gate.
