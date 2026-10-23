# V19 independent dataset-entry audit

## Result

**PASS with a shared-source limitation.** All three fixed-geometry development entries pass the independent input, archive, method, label, and geometry checks listed below. This report does not approve a dataset; A retains that decision.

- Input and archive SHA256 inventories match their pinned manifests and archive checksums.
- GPAW logs show SCF convergence (AgC 59, AgSi 57, AgTi 76 iterations); output coordinates, cell, PBC, atom order, IDs, markers, energies, and finite forces match their registered inputs and summaries.
- Native extrapolated energies match the output and log; free energies are separately recorded and cross-checked.
- Exact/near screen found no hits against train36, six V19 residual input geometries, three public V17 development parents, or the three V18 development frames in `valid.extxyz`.

## Independence boundary and caveats

The three entries use distinct Ag contact configurations and are separated from the audited train/residual and V17/V18 development geometries under the stated periodic comparator. They all derive from the same uploaded Ti3SiC2 COD 9009647 CIF, so they represent one crystallographic source family, not three independent source families. This supports role separation from the checked geometries but does not establish source-diverse generalization or blind test performance. The repository cites COD 9009647 / AMCSD0013252 and Kisi et al. (1998), DOI 10.1016/S0022-3697(98)00226-1; the upstream bytes remain unverified after the recorded HTTP 403.

The fixed PBE/PW 500 eV, 6x6x1, Fermi-Dirac sigma 0.10 eV, xy dipole-corrected labels are method-consistent. Slab k-point, smearing, and vacuum convergence remain unresolved; finite SCF labels do not establish production accuracy. The near-duplicate test only compares equal-composition structures with equal PBC and cells within 1e-5 A, uses 1e-5 A maximum displacement for exact geometry and Ag/framework RMS <=0.15 A for near geometry, and applies no rotation or cell scaling. It is a local screen, not a complete crystallographic symmetry proof.

## Scope and decision state

Only train36 geometry columns, the six residual input geometry files, public V17 development geometry files, and V18 `data/valid.extxyz` geometry columns were screened. The three designated V19 development archives were checked for labels. No test/sealed path was constructed or opened. No dataset was changed or approved, and no preflight, DFT/MPI, MACE inference, or training was run.

`audit.json` contains per-file hashes, per-label evidence, geometry comparisons, and source lineage. `audit.py` reproduces the read-only checks; run it from the repository root with Python plus NumPy, ASE, and SciPy installed.
