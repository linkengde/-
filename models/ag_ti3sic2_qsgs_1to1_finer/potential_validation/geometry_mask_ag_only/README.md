# Ag-only QSGS geometric-mask prototype

This directory contains a finite exploratory implementation of the geometry-only route. The LAMMPS atom system contains Ag only and uses `Ag_u3.eam`; Ti, Si, and C are excluded. The QSGS ceramic occupancy is converted to a signed-distance field and applied as a conservative external penalty force. No Ag–Ti/Si/C pair term or physical interface model is used.

## Model and reproducibility

- The source Ag positions lie exactly on a strained orthorhombic FCC site grid: 73 × 74 × 39 conventional cells, with parent effective lattice constants 4.13638, 4.11082, and 4.09 Å. Preparation preserves the Ag site occupancy and remaps those integer sites to the `Ag_u3.eam` zero-pressure FCC constant, 4.09 Å. The resulting material box is 298.57 × 302.66 × 159.51 Å and keeps all 212,800 Ag atoms. The parent composite's phase ratio is not changed by the coordinate remap, but the Ag-only LAMMPS data does not include bone atoms or bone mass.
- The matching QSGS mask is `potential_validation/relaxed_backbone_candidate/qsgs_morphology.npz` (seed 20261024, 100 × 100 × 50 voxels). The voxel morphology is kept and resized to the remapped Ag box. The wall is open along z and periodic along x/y.
- The bottom Ag atoms initially at z ≤ 20 Å are fixed in the cold-gate runs. The other Ag atoms are mobile. A soft wall clearance of 0 Å and nominal stiffness of 5 eV/Å² were used; stiffnesses 0, 1, 5, and 10 eV/Å² were compared statically.
- The wall energy and force use the same trilinearly interpolated signed-distance field. A central finite-difference check on 123 initially wall-overlapping mobile Ag atoms gave a maximum absolute force error of 5.6 × 10⁻⁹ eV/Å.

Recreate the exported Ag data and distance field with:

```bash
python prepare_ag_only_mask.py
```

The cold test uses the LAMMPS Python interface, an OpenMP EAM build, and `fix external pf/callback`. In the current workspace it was run with:

```bash
LD_LIBRARY_PATH=/workspace/.venvs/mace/lib OMP_NUM_THREADS=4 \
  /workspace/.venvs/mace/bin/python run_ag_mask_cold_gate.py
LD_LIBRARY_PATH=/workspace/.venvs/mace/lib OMP_NUM_THREADS=4 \
  /workspace/.venvs/mace/bin/python run_ag_mask_cold_gate.py --fire
```

The first command is a 2 ps, 300 K NVT screen. The second applies one FIRE minimization capped at 1,000 iterations / 5,000 force evaluations, then performs a 2 ps, 300 K NVT screen. Neither command applies TTM or an external heat source.

## Results and decision

The finite-difference check passed as a numerical implementation check, and the static stiffness sensitivity was recorded. The original-position 2 ps screen retained a 99.974% largest Ag component; the bounded FIRE and subsequent 0.5 ps NVT retained one connected component. These checks support using the mask as an exploratory geometric boundary.

The cold baseline did **not** pass:

- The original-position 2 ps NVT screen fell in potential energy by about 8.2 keV and reached an RMS displacement of 1.30 Å; it was not an equilibrated baseline.
- A bounded FIRE run on the original positions reached the 1,000-iteration cap. It reduced the mobile force 95th percentile to 0.086 eV/Å, but the maximum force remained 1.49 eV/Å. The following 0.5 ps NVT screen ended at 202 K and continued to change the structure.
- The first remapped equilibrium-FCC candidate also did not pass the cold baseline: after 2 ps its potential energy was still falling at about 3.37 keV/ps over the last 0.5 ps, its RMS displacement was 1.26 Å, and its mobile-atom 95th-percentile force was 1.03 eV/Å. That earlier run used a mask interpolation that overcounted the initial surface penetration; the corrected phase-search candidate below supersedes its starting mask count.
- With the nominal soft wall, the direct-extract screens reached minimum signed distance about −0.60 Å. The wall therefore discourages mask penetration but does not impose a hard exclusion.

**No TTM heating or melt-pool run should start from the earlier candidates.** The updated phase-matched candidate and its final bounded cold gate are documented below. Do not interpret this model as a physical Ag–Ti₃SiC₂ interface, adhesion, wetting, interfacial heat transfer, or ceramic thermal response.

## FCC phase search and final bounded cold gate (2026-10-25)

The parent Ag atoms were mapped to the EAM FCC lattice and their rigid in-plane phase was scanned while retaining all 212,800 atoms and the connected network. The best candidate is the original orientation with a small in-plane translation of `(0.06390625, 0.031953125) Å`. Its initial mask overlap improves only slightly, from 174 to 169 atoms (120 mobile); the penetration integral changes from 12.751 to 11.973 Å². The 90°, 180°, and 270° in-plane rotations were rejected because they strongly misalign the FCC sites with the QSGS mask. This bounded scan does not establish a global optimum.

The preparation interpolation at the free z surface was corrected to clamp the continuous z coordinate before interpolation. This removes the previous surface-cell blending error; preparation and force-callback sampling now agree on the candidate's initial 169 mask-overlap atoms.

One finite screen was run on that candidate: FIRE capped at 1,000 iterations / 5,000 force evaluations, followed by 2 ps at 300 K with no heat source or TTM. FIRE hit its iteration cap without reaching the 0.05 eV/Å force tolerance. The Ag network remained one connected component, but the cold gate **failed**:

- The final 0.5 ps potential-energy slope was about `−1,208 eV/ps`.
- The 2 ps RMS displacement from the post-FIRE starting structure was `0.748 Å`.
- Final mobile-atom force p95 was `1.03 eV/Å`; the largest force observed during NVT was `4.86 eV/Å`.
- Atoms inside the soft geometric mask increased from `169` before FIRE to `1,054` after FIRE and `1,088` after NVT.
- The last-window temperature was controlled near 300 K (`300.0 ± 1.7 K`), but this does not offset the continuing energy and force changes.

The wall-force finite-difference error remained small (`3.8 × 10⁻⁹ eV/Å`), so this confirms the force implementation is numerically consistent with its chosen wall energy; it does not validate a physical interface. **Stop this Ag-only geometric-mask route here: do not add more FIRE steps, stiffen the wall, delete additional Ag atoms, or start TTM heating from this candidate.**

## Files

- `prepare_ag_only_mask.py`: exports Ag-only data and creates the signed-distance field.
- `run_ag_mask_cold_gate.py`: finite-difference force check, static wall sensitivity, 2 ps NVT screen, and optional bounded FIRE screen.
- `search_fcc_lattice_phase.py`, `fcc_lattice_phase_search.json`: bounded FCC orientation/translation scan and its results.
- `geometry_mask_metrics.json`: source geometry hashes and preparation statistics.
- `ag_mask_nvt2ps_screening.json`: preserved screen of the original, anisotropically strained Ag coordinates.
- `ag_mask_fire_gate_metrics.json`: bounded FIRE plus 0.5 ps screen of the original, anisotropically strained Ag coordinates.
- `ag_mask_cold_gate_metrics.json`: cold screen of the EAM-equilibrium FCC remapping candidate.
- `ag_mask_phase_fire_gate_metrics.json`: final phase-matched bounded FIRE + 2 ps cold-gate result; the cold gate failed.
- `geometry_wall_callback_history_phase_fire.json`, `log.ag_mask_fire_gate_phase_candidate`: callback trace and LAMMPS log for the final screen.
- `ag_mask_phase_fire_nvt_start.dump.gz`, `ag_mask_phase_fire_nvt_final.dump.gz`: compressed review snapshots from the final screen.
- `Ag_u3_eam_static_reference.json`: zero-pressure bulk reference for the EAM equilibrium lattice constant.
- `ag_mask_cold_start.dump.gz`, `ag_mask_cold_final.dump.gz`: compressed atom snapshots at the start and end of the remapped-lattice cold screen.
- `ag_mask_fire_early_stop_record.json`: earlier FIRE pilot preserved for audit; it stopped early on an overly loose energy criterion and is superseded by the capped run.
- `log.ag_mask_cold_gate_2ps`, `log.ag_mask_fire_gate`: LAMMPS logs for the two final screens.
