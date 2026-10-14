# Fit constraints and identifiability

## Current conclusion

No fit was run. With only the locally auditable data, the TSCBOP parameter groups are not adequately identifiable for a multi-parameter refit. All k-mesh baseline records share one 48-atom geometry; the four displacement structures share that same parent, alter one Ti or C site only along z, and contain no stress. The 114-atom slab result is available only as a report summary in this checkout. The labels therefore do not span independent configurations, strain modes, or surface environments sufficiently to distinguish several coupled Tersoff terms.

## Constraints to preserve in any later pilot

- Keep the raw 27-row file immutable and cite its SHA-256 `729273a6488d75756ed7927432c77c7347a1b115cecbaefae4981bafcec80fd7`.
- Keep all shared triplet values tied exactly as in the source parameterization. Do not split symmetric or shared three-body parameter groups.
- If later authorized, begin with no more than one or two tied parameter groups, in a separate candidate file. All other values stay at source values.
- Publish explicit bounds before optimization. A small ±1% perturbation can be used for finite-difference screening, but is not a fitted bound or evidence of safety.
- Make the objective explicit and dimensionally normalized. A defensible first objective would prioritize atomic force residuals and a separately normalized energy residual; stress should be included only once consistent stress labels, units, and sign conventions are available. Keep native and free-energy conventions distinct.
- Do not fit Ag–Ti/Si/C Morse terms to reduce skeleton residuals.
- Keep the parent bulk geometry family separate from any future validation family. No validation data may be read until its owner releases them.

## Identifiability checks required before fitting

Construct the finite-difference residual Jacobian over the authorized fitting labels and tied parameter groups; scale columns by declared parameter bounds; report singular values, rank under declared numerical tolerance, pairwise column correlations, and leave-one-geometry-family-out changes. A fit should be rejected if several candidate groups produce nearly proportional residual changes or if the Hessian is poorly conditioned. The currently visible labels do not allow these checks to be performed because the engine and direct 114-atom labels are unavailable.
