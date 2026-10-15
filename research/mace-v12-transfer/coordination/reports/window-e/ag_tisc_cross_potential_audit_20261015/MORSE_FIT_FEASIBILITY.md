# calibration-v0 Morse fit feasibility

## Finding

**No parameter fit is authorized by the evidence available here.** A three-pair Morse fit would have up to nine free values (`D0`, `alpha`, `r0` for Ag–Ti, Ag–Si and Ag–C). Current evidence contains one matched full/Ag-only/Ti₃SiC₂-only DFT geometry and six reused-training-family local points without matched phase-only DFT. The 135-atom geometry is a single Si-terminated registry and is not an independent validation family. The ten V19 numerical controls vary k mesh or smearing on one geometry, not physical structures.

## What the matched 135-atom point establishes

It supports one empirical cross-interaction energy anchor and one effective cross-force comparison, with the native and force-consistent energy conventions kept distinct. At this point, the DFT free-energy interaction is −21.378661 eV; the rough LAMMPS run0 decomposition is −47.919378 eV from rounded log output. The archived effective-force RMSE of calibration-v0 against DFT phase-subtracted forces is 1.351321 eV/Å for all atoms and 1.517760 eV/Å for interface atoms within 8 Å. The largest elemental cross residuals are on Si (2.7422 eV/Å) and Ag (2.1496 eV/Å); C is only 0.0946 eV/Å on this one geometry. These are not pairwise independent scores and cannot be taken as the accuracy of the Ag–Si, Ag–Ti and Ag–C functions separately.

One energy and one force geometry cannot determine the curvature and minimum location of all three pair curves. The reused-family local data add only a narrow signed response around their parent structures, and have no matched component DFT needed to strip the known TSCBOP errors. Full-system fitting against those forces would violate the instruction not to let Morse compensate skeleton error. There is no independent family on which to measure generalization.

## Pair-form warning

The existing central pair Morse form cannot satisfy the previously cited 0.05 eV/Å vector-force gate on the two stored Ag–C diagnostic geometries under the fixed intraphase model: even the permissive geometric lower bound that lets each Ag–C neighbor choose an arbitrary independent radial force leaves 0.906–0.909 eV/Å residual. A shared three-parameter Morse curve cannot outperform that lower bound. This does not prove all potential forms fail, nor does it prove the same issue for all Ag–C structures. It does mean that a parameter search over the current form should not be used to promise that gate. C must decide whether the Ag–C functional form or acceptance metric needs revision before a calibration attempt.

## Identifiability status

No cross-parameter Jacobian/SVD or multi-geometry independent sensitivity report is present in the reviewed cross archives. It would be incorrect to infer parameter rank from the species-level force RMSE table or from a single geometry. The six signed points do not provide six independent families, and no exact DFT phase-resolved cross-force targets are available for them. Therefore:

- `parameter_identifiability = NOT_DEMONSTRATED` for nine independent Morse parameters;
- `cross-force targets = PARTIAL` (one geometry has matched phase differences);
- `energy-distance coverage = INSUFFICIENT` (one matched interaction-energy point, plus reused local perturbations not equivalent to a normal scan);
- `independent validation = UNAVAILABLE` in the reviewed data;
- `trial fit = NOT RUN`.

A future trial should first freeze one parameterization strategy. At minimum compare the nine-parameter model against a constrained version with alpha frozen to the existing common 2 Å⁻¹ value and only pair-specific depth/equilibrium distance varied (six parameters). Do not assume that freezing alpha makes the model identifiable; compute the force-and-energy response Jacobian on the frozen training structures and require full rank with acceptable conditioning before optimizing. The Ag–C lower-bound warning remains in force.

## Minimum conditions before a trial

1. C freezes a specific set of physically admissible interface structures, target DFT settings, role assignments and the Ag–C force/functional-form acceptance decision.
2. TRAIN includes a fixed-geometry energy/force curve with at least three normal separations (repulsive, near minimum, and attractive side) at two non-equivalent lateral registries. The structures must expose meaningful contact-distance variation for every pair being fitted; if Ag–C remains far from contact on the chosen termination, include a physically justified C-accessible termination or do not claim Ag–C calibration.
3. A distinct parent surface/registry family is held out from all fitting and parameter selection. It must cover the target pair environments. It cannot be a different k mesh, smearing, or perturbation fragment of a training parent.
4. All training targets have SCF-converged PBE/PAW energy/free energy and forces with one frozen numerical method. Stress is retained where generated, never fabricated. For cross-only targets, calculate full and matched separated phases on the same cell/coordinates. Reuse existing phase references only after exact component geometry, cell, PBC, atomic-ID mapping and rigid-transform invariance are verified.
5. Evaluate analytic or finite-difference sensitivities and rank/conditioning before fit. Fit DFT interaction energy/free energy differences and phase-subtracted forces with consistent weights; do not fit absolute energy offsets across unlike compositions. Keep TSCBOP and Ag EAM fixed.

Even after a successful trial, this does not establish production readiness. Small-model static, controlled minimization and 300 K no-heat interface checks would need a separately versioned candidate after C approves it. No such candidate exists from this task.
