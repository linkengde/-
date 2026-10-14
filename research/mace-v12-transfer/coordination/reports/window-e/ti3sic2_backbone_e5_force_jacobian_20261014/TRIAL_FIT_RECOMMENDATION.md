# Stage E5 trial-fit gate

**Recommendation: START_TRIAL_FIT=NO.** This stage completed the authorized static force-response analysis only. The Jacobian shows two unexcited parameters, high coupling among active fields, and no independent DEV/TEST family. The current 48-atom perturbations all belong to one parent; the 114-atom TSC-only slab is derived and is not an independent clean-surface validation set. No roles were changed.

The smallest numerically separated two-column candidate is `A_Ti-C, B_Ti-Si` (condition 3.288 after unit-column scaling), but the response matrix is not a substitute for a balanced DFT objective or held-out validation. Before any separately authorized fit, obtain Si–C-active geometries, independent volume/shear/Si-displacement and clean-surface data as E3 specified, then freeze TRAIN/DEV/TEST by source family. If a fit is later approved, hold all other Tersoff blocks fixed and regularize toward the original parameterization; the ±1% probe itself is not a bound.

No TSCBOP, EAM, Morse, DFT labels, or frozen structures were modified. No fitting/optimization, DFT, MD, TTM, or MACE run was started. Candidate parameters remain diagnostic only and are not production-ready.
