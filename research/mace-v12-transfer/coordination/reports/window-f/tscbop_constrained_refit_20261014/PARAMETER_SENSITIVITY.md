# TSCBOP parameter sensitivity

## Status

`NOT_TESTED`. No new force evaluations were possible because ASE and LAMMPS are absent in this task container. The audit did not install them, invoke GPAW, run MD, or edit the raw parameter file. Consequently there are no numerical finite-difference sensitivities or parameter-to-metric effect estimates in this report.

## Parameter families that can affect the solid framework

The 27-line file has the standard Tersoff triplet fields after its three element labels: `m, gamma, lambda3, c, d, h, n, beta, lambda2, B, R, D, lambda1, A`. In the Tersoff energy, `A` and `lambda1` shape the repulsive radial term; `B` and `lambda2` shape the attractive radial term; `beta` and `n` control bond-order response; `gamma`, `lambda3`, `c`, `d`, and `h` control coordination and angular response; `R,D` set the cutoff envelope. All can influence Ti/C/Si forces, but through overlapping terms rather than one-to-one element-force channels.

The current file has only one exact value each for `m`, `n`, and `beta`; six for `gamma`, `c`, `d`, `h`, `lambda2`, `B`, `R`, `lambda1`, and `A`; seventeen for `lambda3`; and five for `D`. These column counts are descriptive only. They do not establish independent physical parameters: many triplet entries share values and are coupled by the directed triplet table. Any future perturbation must operate on explicitly tied groups and preserve the published table symmetries. A row-by-row unconstrained adjustment would violate the requested shared three-body structure.

## Safe screening design if the audit gate is later cleared

1. Freeze and hash the raw file; create a separate candidate copy.
2. Validate the 27-row parser and the three-element mapping against a known original-file run0 before perturbing parameters.
3. Define tied parameter groups from the source table, not from per-row convenience. First use one-group-at-a-time small perturbations (for example ±1% as a screening step, not a fitted bound) and report changes in force RMSE, native/free energy separately, stress if available, and the signed Ti/C displacement response.
4. Only after screening, choose at most one or two weakly correlated groups for a bounded pilot. Do not vary all radial, angular, cutoff, or bond-order terms together. Keep every nonselected value at the raw-file value.
5. Compute a finite-difference Jacobian of the full labeled residual vector and its singular values/correlation matrix. Reject a multi-parameter proposal if the selected columns are effectively collinear at the available data precision.
6. Do not tune Ag cross terms. Do not use the same-parent displacement pair as independent validation.

This is a proposed protocol, not an executed sensitivity analysis. No candidate values or metric deltas were generated.
