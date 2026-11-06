# Ag70 small-pool visualization at 4.52 ps

- Source: existing Ag70, 1×, Case A checkpoint `caseA_ag70_1x_final_4p52ps.npz`; no new heat simulation was run.
- Time: `4.520000 ps`; peak temperature `1264.383 K`.
- Main connected cluster uses the solver's `f_liquid >= 0.5` threshold: `65` native voxels.
- Threshold-mask geometric volume: `1.864973 nm³`; liquid-volume integral within that cluster (`sum(f_liquid) × voxel volume`): `1.838613 nm³`.
- Axis-aligned voxel extent: `2.397080 × 1.811145 × 0.951669 nm`.
- There is one `f_l >= 0.5` connected cluster at this snapshot. The 4×4×4 global context is only for location; the right inset and reported metrics use native voxels.
- This enthalpy model contains no arc-root plasma, material evaporation, liquid flow, or splashing.
