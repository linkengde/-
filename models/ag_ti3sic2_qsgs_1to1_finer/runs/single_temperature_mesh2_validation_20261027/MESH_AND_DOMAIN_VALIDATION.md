# Mesh and domain validation status: Ag 70:30 vol%

The selected geometry is the existing QSGS growth-rank realization thresholded to Ag 69.9996 vol% / Ti₃SiC₂ 30.0004 vol%. It was not re-randomized. The thermal domain is 100×100×50 native voxels (approximately 3 Å spacing), periodic in x/y and open/insulated in z.

## Existing target-composition check through 0.2 ps

The 1× and 2× cases use the same binary geometry, exact Gaussian-plus-depth source integration, Case A, initial 300 K, and unchanged source/material inputs. At 0.2 ps their peak temperatures are 748.6093 and 723.1792 K, respectively, and both have zero liquid fraction. Exact source powers agree to roundoff; relative total energy residuals are approximately −1.02×10⁻¹⁴ and 4.94×10⁻¹³. Block-averaged 2× versus 1× Ag temperatures have MAE 0.404 K and RMSE 1.415 K, but the peak differs by 25.43 K and moves closer to the source centre. Thus the coarse local maximum is mesh-sensitive; the no-melt conclusion at 0.2 ps agrees. This is not a melt-pool convergence result.

## New 1× melt-window check through 1.0 ps

The 1× and 2× 0.2 ps checkpoints were restored and their metadata/field hashes and energy ledgers accepted by `load_checkpoint()`. A 1× continuation from 0.2 to 1.0 ps took 904.39 s. Peak temperature increased from 748.61 K at 0.2 ps to 906.21 K at 1.0 ps; every sampled state remained fully solid. Relative energy residual at 1.0 ps was 5.53×10⁻¹³ and the final checkpoint round-tripped. Since no liquid formed, neither pool-size convergence nor cooling-time convergence can be tested in this window.

The 2× grid is currently being continued from 0.2 to 1.0 ps. Its same-time comparison will be inserted after completion. The prior 0.2 ps comparison showed 2× peak temperature 25.43 K below 1×, while Ag-weighted block-averaged MAE/RMSE were 0.404/1.415 K. That early-time result motivates the matched-time check; it does not establish the 1 ps behavior.

## Domain-size validation

The previous 3×3 tile study belongs to the old 41.6654 vol% Ag geometry and cannot validate this 70:30 composition. No 70:30 3×3 heating result is claimed. Before any such run, estimate memory and runtime, tile the mask only in x/y, keep exactly one central source, verify exact integrated source power and energy balance, and compare central-cell-region fields rather than total nine-cell liquid volume. If time is limited, leave it as a costed proposal.

## What remains necessary

- Compare 1×/2× at shared times around first liquid and first connected multi-voxel liquid; use physical volume for clusters and report field differences after block averaging.
- A single pair at one voxel scale interval is a mesh-sensitivity check, not Richardson convergence. A further refinement or local error study may be needed if onset or pool dimensions differ materially.
- Check timestep sensitivity at the phase transition and in source-off cooling.
- Domain periodicity has not been checked for the target composition. A 3×3 repeat tests only lateral-image interactions for this same realization, not independent-QSGS morphology variability.
- The approximately 3 Å native voxel scale limits minimum-neck and interfacial morphology interpretation. Grid-scale segmentation can change connected-liquid onset.
