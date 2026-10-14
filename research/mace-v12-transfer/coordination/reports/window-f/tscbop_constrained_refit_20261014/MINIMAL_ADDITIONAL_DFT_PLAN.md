# Minimal additional DFT evidence plan (proposal only)

No DFT was run or launched for this task. The following files already exist as unexecuted proposals under `research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_pure_phase_controls/inputs/`; they are not labels:

- `Ti3SiC2_volume_m02_PROPOSAL.extxyz` and `Ti3SiC2_volume_p02_PROPOSAL.extxyz`;
- `Ti3SiC2_shear_m005_PROPOSAL.extxyz` and `Ti3SiC2_shear_p005_PROPOSAL.extxyz`.

If separately authorized, execute these four fixed-geometry PBE/PAW points with the same 48-atom parent, 500 eV cutoff, 0.10 eV smearing, and converged mesh. Archive native and free energies separately, forces, stress with explicit units/sign/Voigt convention, convergence evidence, and exact input/output/setup hashes. These would constrain a first bulk volume/shear response; they would not independently validate the potential or establish finite-temperature stability.

Before a fit, also restore the raw 114-atom TSC-only archive and its source manifest so the interface-derived diagnostic can be audited directly. For independent validation, obtain a separately sourced bulk or distinct-displacement parent and freeze its geometry and labels outside the fitting workflow. A clean surface or interface-neighbor label should be added only if that environment is part of the validation scope, with its provenance and geometry role reviewed in advance.

This is a minimal evidence recommendation, not authorization to execute DFT.
