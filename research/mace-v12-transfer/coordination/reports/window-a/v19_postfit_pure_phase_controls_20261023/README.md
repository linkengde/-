# V19 post-fit pure-phase diagnostic

## Finding

The V19 epoch-79 model reproduces the selected Ag fcc control, but retains a large Ti₃SiC₂ host-force error on the existing pure-phase controls. This confirms that the V19 failure is not limited to Ag–Si cross interactions. It does **not** turn these controls into independent validation or training data.

All 26 calculation archives in `pbe_interface_v19_pure_phase_controls/calculations/` were checked in this read-only pass: SCF convergence, per-archive SHA256 entries, and recorded verifier checks all pass. The selected seven inputs/results were also checked for exact atom order, positions, cell, PBC, and source hashes before inference.

| Existing numerical control | Force-vector RMSE (eV/Å) | Ti RMSE | C RMSE | Si RMSE |
|---|---:|---:|---:|---:|
| Ag fcc, 8×8×8 | 0.0000173 | — | — | — |
| Ti₃SiC₂ baseline, 10×10×2 | 3.00887 | 3.58702 | 2.80352 | 0.00135 |
| Ti₃SiC₂ baseline, 10×10×4 | 3.00829 | 3.58637 | 2.80291 | 0.00076 |
| Ti₃SiC₂ Ti4f ID3, −0.20 Å | 3.01609 | 3.59839 | 2.80494 | 0.01039 |
| Ti₃SiC₂ Ti4f ID3, +0.20 Å | 3.00041 | 3.57491 | 2.79952 | 0.01217 |
| Ti₃SiC₂ C4f ID9, −0.20 Å | 3.00622 | 3.58482 | 2.79924 | 0.00331 |
| Ti₃SiC₂ C4f ID9, +0.20 Å | 3.01071 | 3.58876 | 2.80614 | 0.00314 |

These force-vector errors are model-minus-DFT residuals, not raw predicted force magnitudes. The Ti₃SiC₂ reference structures have small nonzero forces; the comparison uses the registered DFT force arrays. This is a diagnostic on numerical-control geometries and must not be presented as an independent model score.

## Reuse and source boundary

The six Ti₃SiC₂ results are valuable for diagnosing backbone response, but every calculation record explicitly assigns `numerical_pure_phase_control_not_training_or_independent_validation`. They derive from COD 9009647, the same parent source family as V19 dev3. Do not add them to training while retaining dev3 as a source-independent check; that would leak the current development source family into fitting.

The repository contains two byte-different CIF files for Ti₃SiC₂, but both identify COD 9009647, the same 1998 Kisi et al. source, and the same crystallographic entry. Their byte hashes differ because one is the user's original copy and the other is a reloaded COD copy; this does not provide an independent source family. No second qualified Ti₃SiC₂ CIF parent is present in the project checkout. A public Materials Project OPTIMADE query was attempted, but the network tunnel returned HTTP 403.

The Ag pure-phase control passes, and the existing Ag–Ti₃SiC₂ interface errors are dominated by Ti/C forces. New data should therefore prioritize a source-verified second Ti₃SiC₂ parent and force-rich host/backbone configurations, while keeping the present COD 9009647 dev source out of fitting. Before any new fit, reconcile train36's Gamma-point, fully periodic labels with the slab dev method and audit inherited energy keys.

No role was changed, no DFT/MPI or training was started, and no test/sealed data were opened.

## Files

- `metrics.json`: archive/input hashes and per-atom inference results.
- `infer_v19_pure_phase_controls.py`: reproducible read-only archive audit and inference.
- `SHA256SUMS.txt`: hashes for this bundle and its referenced selected data/model artifacts.
