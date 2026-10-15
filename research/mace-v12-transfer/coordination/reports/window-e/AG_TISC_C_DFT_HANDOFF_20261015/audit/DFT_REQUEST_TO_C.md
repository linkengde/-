# DFT request to C: minimum cross-interaction evidence

Status: **request for C review and structure/role freeze only**. E did not launch DFT or assign TRAIN/DEV/TEST roles.

## Why more data are needed

The reviewed package has one matched full/Ag-only/Ti₃SiC₂-only fixed-geometry DFT point (the 135-atom Si-terminated interface), plus six converged V19 local references from already-used training families. The latter do not have matched phase-only DFT and cannot safely be used to fit Ag–X cross forces with the known TSCBOP force error held fixed. Existing points do not form a normal-distance curve, independent registry set, or untouched validation family. The C bundle's ten extra cases are numerical controls of one geometry.

The archived Ag–C geometric force lower bound also shows that the current central Morse basis cannot pass the existing 0.05 eV/Å vector gate on two stored Ag–C inputs, even if each pair force amplitude is freed independently. Before authorizing any Ag–C fit, C should decide whether the current Morse form and that gate remain the project requirement.

## Proposed minimum staged acquisition

**Freeze before acquisition:** select exact input structures, periodic cell, terminations, registries, contact-distance targets and a single DFT configuration. Preserve all source roles and hashes. Keep the existing 135-atom structure/results immutable.

**TRAIN coverage floor:** six fixed-geometry full-interface DFT single points arranged as two non-equivalent lateral registries × three normal separations (compressed/repulsive side, near the local minimum, and expanded/attractive side). The existing 135-atom full-interface result may count as one point only if C places its measured gap and registry on this frozen grid; otherwise it remains a separate training-family anchor and all six points are new. Include each target Ag–Ti, Ag–Si and Ag–C environment at a non-negligible distance. If the selected physical termination never activates one pair sufficiently, add a pair-relevant, justified termination series rather than trying to infer that pair from a distant contact.

For the current 135-atom parent and rigid translations/separation changes that leave each isolated phase internally identical, its existing Ag-only and TSC-only calculations may be reused only after exact cell/PBC, phase coordinates modulo the specified rigid transform, and source-atom mapping are verified. For any changed slab/termination/strain where this invariance does not hold, add a matching neutral Ag-only and TSC-only reference at the same cell and frozen component geometry. For every target geometry, retain full-interface force-consistent free energy, native energy, forces, stress (if calculated), SCF result, cell/PBC and input/result hashes. Use `E_interaction = E_interface − E_Ag − E_TSC`; keep native and free-energy differences distinct.

**Validation family:** three fixed-geometry full-interface single points at three corresponding gaps for one independently constructed parent surface/registry family not used for training or parameter selection. It must test the pair environments intended for use. Acquire the two matched phase-only references for that validation family unless exact, verified phase invariance permits reuse. C owns the final role and test-set seal.

This is a **floor of nine full-interface geometries** across train and held-out validation. Under the narrow assumption that the existing 135-atom result occupies one TRAIN grid slot and its two phase references are valid for the training scans, the floor implies five new TRAIN full-interface points, three held-out full-interface points, and up to two phase-only references for the new validation family: **at least ten new SCF single points total**. If the 135 geometry does not occupy a grid point, or if the training/validation phase structures change, add those missing full or component jobs. If Ag–C (or another target pair) lacks meaningful distance coverage, add a three-point pair-relevant termination curve and matched components; the ten-job floor is then insufficient.

This lower bound is not a promise that the resulting nine Morse parameters will be identifiable. After the frozen geometries are assembled, first run only the classical static sensitivity/geometry-rank check on temporary copies and check the analytic force/energy Jacobian rank/conditioning. If a pair remains unidentifiable, C should approve only the smallest additional DFT set needed for that pair. Do not duplicate the existing six signed references or ten k/smearing controls.

## Frozen calculation and data requirements

- Fixed atoms and cell; no relaxation, MD, or thermal trajectory.
- PBE/PAW, GPAW 26.7.0, 500 eV, Γ 1×1×1, FD 0.10 eV and the archived 1e-5 SCF thresholds as the initial common-method proposal; C may issue a corrected locked method before authorization.
- Store native total energy and force-consistent free energy separately, full per-atom three-dimensional forces, stress if requested, converged SCF steps/status, atom IDs/order, cell/PBC, element counts and SHA256 manifest.
- Component-only calculations must be neutral, same cell/PBC/coordinates for the retained atoms, and share the full-system functional, PAW setup and numerical settings. Record valence electron count and occupations; flag rather than conceal any spin/odd-electron concern.
- Do not use absolute DFT total energies across compositions. Use matched interaction-energy differences and matched phase-subtracted force differences. Do not make a Morse term cancel TSCBOP's intraphase error.

## Requested C decisions

1. Approve/freeze the target interface family, allowed terminations, registry shifts and gap grid.
2. Decide whether the current pairwise Morse form and 0.05 eV/Å force gate remain acceptable in light of the Ag–C lower-bound evidence.
3. Freeze roles: proposed TRAIN = six-point two-registry gap set; proposed validation = separate three-point family. Do not alter existing dataset roles without registry review.
4. If the acquisition floor is accepted, issue a new DFT authorization tied to exact inputs/configuration SHA. Until then, E remains HOLD and will not execute calculations or fitting.
