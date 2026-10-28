# Ag 70 vol% / Ti₃SiC₂ 30 vol%: 1× heating and cooling result

**Status:** the requested 1× Case A heat–melt–source-off cooling branch completed. A short-lived connected liquid cluster formed and returned to the all-solid state. This is a transient numerical pool under the stated idealized model, not evidence of a stable weld pool or an actual Ag–Ti₃SiC₂ weld. The matched 2× melt-window comparison remains pending.

## Run and verification

The calculation continued from the existing 1.00 ps checkpoint with the original source, properties, geometry, and Case A boundary conditions. The environment disconnected after a valid 2.20 ps checkpoint; that checkpoint was read back and the same run resumed from it. The heating log covers the full 1.00–4.52 ps continuation. No Case B or separate double-temperature task was touched.

The continuation used refinement 1, a 5.0 ps maximum time, a 2300 K numerical temperature guard, and 0.30 ps of follow-up after the first multi-voxel liquid cluster. It stopped at 4.52 ps because the follow-up completed. The source was then turned off and the saved 4.52 ps state was continued with a 0.5 ps cooling limit, 0.02 ps segments, and 0.002 ps samples. Cooling stopped at 4.56 ps after the all-solid condition was sampled.

The final heat-on checkpoint passed the continuation script's restore-only validation: time 4.52 ps, peak 1264.38 K, Ag liquid fraction 1.8923×10⁻⁴, 49 tracker history rows, and relative energy residual −1.61×10⁻¹². The terminal cooling checkpoint arrays were read back successfully; the temperature, enthalpy, and liquid fields are finite, and its liquid fraction is zero.

## Melt and cooling metrics

| Event or measure | Result |
| --- | ---: |
| First nonzero Ag liquid fraction | 4.12 ps |
| First connected component | 4.19 ps; one voxel |
| First multi-voxel liquid cluster | 4.22 ps; two voxels, 57.38 Å³ |
| Heat source switched off | 4.52 ps |
| Peak temperature at switch-off | 1264.38 K |
| Ag liquid fraction at switch-off | 1.8923×10⁻⁴ (fraction of Ag volume) |
| Largest connected cluster at switch-off | 65 voxels, 1864.97 Å³ (1.865 nm³) |
| Cluster centroid at switch-off | (153.75, 147.47, 154.83) Å |
| Cluster depth from top surface | 3.78 Å |
| Lateral span | 23.97 Å × 18.11 Å |
| Main track’s sampled continuous duration | 0.31 ps |
| Last sampled nonzero liquid state after source-off | 4.546 ps |
| First sampled all-solid state after source-off | 4.548 ps |
| Sampled source-off-to-all-solid interval | 0.028 ps |

The largest track persisted through the 0.30 ps follow-up. Its volume increased from a few voxels to 65 voxels, and the tracked liquid components later split and merged. The 0.31 ps tracking duration is a sampled connectivity measure; it is not a general stability threshold. At switch-off, the model liquid volume was about 0.0189% of the Ag volume, or about 0.01325% of the whole composite volume at 70% Ag.

The reported solidification interval is bounded by samples: liquid remained at 4.546 ps and was absent at 4.548 ps. It should be read as a roughly 26–28 fs interval after switch-off at the current 2 fs sample spacing. The insulating domain does not lose energy to its surroundings; cooling here is internal redistribution and latent-to-sensible energy conversion, not an environmental heat-loss time.

## Energy balance

At switch-off, cumulative source energy was 7.012142437008×10⁻¹⁵ J and total enthalpy was 7.012142436997×10⁻¹⁵ J. The recorded residual was −1.1320×10⁻²⁶ J (relative −1.614×10⁻¹²). Total enthalpy remained unchanged through the source-off branch, within the same residual. The run did not approach the 2300 K numerical guard; its maximum was 1264.38 K.

## Resolution check and limits

The available 2× target-composition checkpoint ends at 0.35 ps and is fully solid. The earlier 1×/2× comparison at 0.35 ps has peak temperatures 809.70/777.19 K, a 32.51 K difference; it cannot validate a molten-pool metric. Reaching the 1× multi-voxel pool time from the saved 2× state would require about 3.87 ps more simulation. Using the previously measured 2× pace of roughly 97 minutes per 0.2 ps gives an extrapolated cost of about 31 hours, so no low-value partial 2× run was launched. The 1× melt-pool size and lifetime are therefore not mesh-converged.

The model still uses candidate thermophysical properties not fully source-verified, a single temperature field, a prescribed Gaussian-plus-depth-exponential source (30 Å radius, 20 Å absorption depth, about 1.551 mW absorbed power), periodic lateral boundaries, and Case A’s zero-conductivity/zero-heat-capacity Ti₃SiC₂ skeleton with insulated exterior. Vaporization, finite Ag/Ti₃SiC₂ interface conductance, scaffold heat storage/conduction, and real contact-current coupling are absent. These results are exploratory model responses and cannot establish real contact welding or comparative anti-welding performance.

## Artifacts

- Combined heating/cooling metrics: [MELT_FREEZE_METRICS.csv](MELT_FREEZE_METRICS.csv)
- Temperature, liquid fraction, and largest-cluster history: [ag70_melt_freeze_evolution.png](ag70_melt_freeze_evolution.png)
- Spatial snapshots: [first multi-voxel cluster, 4.22 ps](ag70_pool_first_multivoxel_4p22ps.png), [heat-off, 4.52 ps](ag70_at_heat_off_4p52ps.png), [cooling, 4.54 ps](ag70_cooling_4p54ps.png), [all solid, 4.56 ps](ag70_after_solidification_4p56ps.png)
- Reproducible summary/plot generator: [summarize_ag70_melt_freeze.py](../summarize_ag70_melt_freeze.py)
- Full heating segments and cluster tracks: [segments.jsonl](../ag70vol_resume_1x_from1ps_20261028/segments.jsonl) and [heating.clusters.csv](../ag70vol_resume_1x_from1ps_20261028/heating.clusters.csv)
- Cooling history and cluster tracks: [cooling.csv](../ag70vol_cooling_20261028/cooling.csv) and [cooling.clusters.csv](../ag70vol_cooling_20261028/cooling.clusters.csv)
- Final heat-off and cooling restart states: [4.52 ps heat-off state](../ag70vol_resume_1x_from1ps_20261028/checkpoints/caseA_ag70_1x_final_4p52ps.npz) and [4.56 ps all-solid state](../ag70vol_cooling_20261028/checkpoints/caseA_ag70_1x_cool_final_4p56ps.npz)

The run directories retain the intermediate restart checkpoints locally. The 2× matched melt-window check is still needed before making a resolution-independence claim.
