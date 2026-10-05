# v15 force and separation-force correction plan

## Evidence

V14 trained to optimizer epoch 79 and selected epoch 71. Its training configuration already used forces_weight=1000 and energy_weight=100; force weight was ten times the energy weight. Therefore, changing only the force weight is not the evidence-led first fix.

The provisional gates are energy MAE <= 10 meV/atom, force-vector RMSE <= 0.05 eV/Angstrom, and marked-pair separating-force absolute error <= 0.10 eV/Angstrom. On the three V14 blind holdouts, energy MAE passed (Ag-C 3.23, Ag-Si 4.62, Ag-Ti 1.68 meV/atom), while force RMSE failed for all (0.154, 0.118, 0.147 eV/Angstrom). Separating-force error passed for Ag-C (0.061) and failed for Ag-Si (0.235) and Ag-Ti (0.193 eV/Angstrom).

The gap distances are already represented closely in V14 training:
- Ag-C holdout 2.249 Angstrom; V14 training points 2.224 and 2.260 Angstrom.
- Ag-Si holdout 2.439 Angstrom; V14 training points 2.440 and 2.463 Angstrom.
- Ag-Ti holdout 2.564 Angstrom; V14 training points 2.588 and 2.600 Angstrom.

So the large Ag-Si/Ag-Ti separation-force errors cannot be attributed to missing only the scalar pair distance. Local registry, angular coordination and the surrounding force response are likely contributors.

The V14 force-localization report shows:
- Ag-C: local RMSE 0.176 versus 0.094 eV/Angstrom outside 4.5 Angstrom; marked-pair separation error passes.
- Ag-Ti: local RMSE 0.178 versus 0.083 outside; maximum single-atom error 0.465 eV/Angstrom at C#480805, near Ti#479027 and Ag.
- Ag-Si: local RMSE 0.127 and outside RMSE 0.108; error is not confined to the marked contact shell.

Window B's geometry audit found each current v15 candidate is only a small perturbation of a frame already in V14 training (same-ID positional RMS 0.029, 0.052 and 0.044 Angstrom for Ag-C, Ag-Si and Ag-Ti). These labels can be useful training additions but cannot serve as independent v15 blind tests. No unused distinct registry input is already present in the repository.

All current geometries are roughly 48-51 atom cluster motifs in approximately 28 Angstrom cubic cells. This evidence does not establish accuracy for an extended interface, hot structure or liquid Ag configuration.

## Dataset-integrity finding from the V14 audit (2026-10-05)

The V14 frozen AgTi test geometry is an exact train/test geometry duplicate. `AgTi_2p7A_LCAO` in the V14 training set and `AgTi_2p70_periodic_PW_PBE_energy_force_holdout` in the test set have the same 49 atom IDs, elements and coordinates (same-ID positional RMS and maximum are both 0). The training frame is explicitly marked force-only LCAO, while the test frame is a PW-PBE reference. Comparing the two force arrays at this identical geometry gives a vector RMSE of 0.0852 eV/Angstrom and a largest per-atom difference of 0.2013 eV/Angstrom. This test therefore is not an independent geometry test, and the two reference-label sets disagree beyond the provisional 0.05 eV/Angstrom force gate. The cause cannot be assigned to basis, periodic-boundary, or convergence settings separately without the original matched calculation records.

V14 has four legacy force-only LCAO training frames in total (two AgC and two AgTi); it has 31 training structures but only 27 energy labels. Treat these as a reference-consistency risk until their source settings and matched PW-PBE labels are audited. Do not use the duplicate AgTi frame as independent evidence of generalization. Keep all V14 files and scores unchanged as historical evidence; apply the corrections to v15 data and split design.

The V9-V15 history report added at commit `c48ae13` conflates two AgTi structures in its independence note. `AgTi_registry_holdout_v14_01` is **not** an exact geometry duplicate of `AgTi_2p7A_LCAO`: it is derived from the frozen 2.70 Angstrom parent but uses a 2.60 Angstrom target contact; its nearest V14 training geometry has same-ID positional RMS 0.4153 Angstrom and maximum displacement 0.5746 Angstrom. The exact duplicate is instead the separate frozen-test frame `AgTi_2p70_periodic_PW_PBE_energy_force_holdout`, which exactly matches `AgTi_2p7A_LCAO`. The reported historical scores on the V14 registry holdout remain scores on a distinct geometry; its broader independence is still limited by shared parent-motif lineage and by unverified V9-V11 training membership.

## Changes for v15

1. Finish and archive the existing three residual-targeted DFT labels. If verification passes, add them to v15 training only. Do not count them as validation or blind tests.
2. Audit v15 train/valid/test and holdout geometries for exact and near duplicates by persistent atom ID, element, and coordinates before training. The V14 AgTi frozen test is an exact duplicate of a training geometry and must not be counted as an independent test.
3. Resolve reference consistency for the four legacy LCAO force-only training frames. Use matched PW-PBE labels where available; otherwise exclude them from the PW-PBE v15 objective until their method/boundary/convergence provenance is established. Never keep a geometry in a blind split when another reference label for that same geometry is in training.
4. Prepare new label-free holdout geometries at a genuinely different Ag registry for Ag-C, Ag-Si and Ag-Ti. Use controlled rigid Ag translation plus small, recorded local displacements; reject candidates with unreasonable short contacts or near-duplicate geometry. Freeze and hash these inputs before DFT scoring. These test new local registries, not independent material morphologies.
5. Keep each new holdout out of v15 training, validation, checkpoint selection and hyperparameter selection. Preserve V14's frozen tests and previously scored holdouts as historical regression evidence, with the AgTi exact-duplicate limitation clearly annotated.
6. To reduce all-atom force error, acquire training labels for multiple local environments around the Ag-C and Ag-Ti hot shells, and include broader neighboring/framework distortions for Ag-Si because its error extends beyond the contact shell.
7. To reduce separating-force error, evaluate the marked-pair force component at new registries. Add a short pair-distance scan around a new registry only if the independent screen still shows separation-force error; the existing distance points are already close to the V14 test distances.
8. In the first v15 comparison, keep the MACE architecture, foundation model, seed, optimizer and force/energy weights fixed at the V14 values. This isolates the effect of new data. Do not simply increase force weight or epoch count without validation-curve evidence. If forces remain poor after broader data coverage, run a controlled weight ablation on a validation set while leaving the fresh blind test untouched.
9. Report, per interface and per structure: energy error, all-atom force-vector RMSE, local and outside-shell RMSE, largest atom residual with neighbors, and marked-pair separating-force error. A passing single structure per interface remains a screening result, not broad transfer validation.

## Current queue and gates

AgC_residual_shell_v15_01 converged in 60 SCF iterations and passed archive verification. AgSi_residual_shell_v15_01 converged in 49 iterations and passed archive verification; it was published with its compact log and manifest at commit `9c5c670`. AgTi_residual_shell_v15_01 remains queued. Do not restart completed jobs.

B's label-free registry screen nominated AgSi and AgTi proposals for A review. Both pass its geometry/duplicate checks, but they are translations of small motifs already present in v14 training, so any DFT labels should be described as new local-registry holdouts only. The AgC proposal was withheld because its unchanged framework C#604587-Ti#604590 pair is 1.7413 Å. Before rejecting it, assess that pair against species-specific Ti-C framework distances and connectivity; a universal all-pair cutoff alone does not establish that a bonded framework pair is unphysical. A fresh AgC proposal must still pass a chemically appropriate geometry screen and the duplicate audit.

The next gates are: finish/verify the acquisition queue; generate and freeze/hash new registry holdouts; build and audit the v15 dataset; train once with the controlled baseline; score old regression and fresh blind structures; decide PASS, FAIL or UNDETERMINED. No production long MD or TTM before the force and separation-force screens pass.
