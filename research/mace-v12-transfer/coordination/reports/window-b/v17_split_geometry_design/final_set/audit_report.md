# V17 validation/test geometry design (label-free)

## Current decision

Window B prepared three development-validation proposals and six withheld-test geometry proposals, all before V17 training. Every structure is geometry-only and preserves atom count, atom order, persistent IDs, cell and PBC. The role manifest is proposed as frozen; A must review and accept it before training. Withheld-test labels must not be used for checkpoint or parameter selection.

## Candidate screen

| Candidate | Proposed role | Atoms | Marked pair (ID; distance Å) | Exact/near inventory hits | Species-pair checks | File SHA256 |
|---|---|---:|---|---:|---:|---|
| AgC_v17_dev_validation_01 | development_validation | 48 | 688433–604804 (2.6597) | 0 | 10 | f0d8bb9d0652… |
| AgC_v17_test_transverse_framework_01 | withheld_test | 48 | 688433–604804 (2.8023) | 0 | 10 | e58bc3e66618… |
| AgC_v17_test_distance_01 | withheld_test | 48 | 688433–604804 (2.0600) | 0 | 10 | 5e3d40c3c79f… |
| AgSi_v17_dev_validation_01 | development_validation | 51 | 683609–366934 (2.8877) | 0 | 10 | 45333ff0759b… |
| AgSi_v17_test_transverse_framework_01 | withheld_test | 51 | 683609–366934 (2.9198) | 0 | 10 | 74ed395968a9… |
| AgSi_v17_test_distance_01 | withheld_test | 51 | 683609–366934 (2.2651) | 0 | 10 | 78af998abb01… |
| AgTi_v17_dev_validation_01 | development_validation | 49 | 687247–479027 (3.2230) | 0 | 9 | 7c83fad1f8b7… |
| AgTi_v17_test_transverse_framework_01 | withheld_test | 49 | 687247–479027 (2.9727) | 0 | 9 | 06a0c3ea33d8… |
| AgTi_v17_test_distance_01 | withheld_test | 49 | 687247–479027 (2.7630) | 0 | 9 | 247cc873f175… |

The comparison inventory contains 197 frames from 18 prior model-split files and 38 geometry frames from 38 published DFT input files. This includes all 35 V16 training rows and the eight hash-verified planned V17 training inputs. Historical probe inputs are included via their input paths only.

Exact duplicates use persistent IDs when available, equal PBC/cell and maximum minimum-image displacement ≤1e-5 Å. Near duplicates use Ag registry RMS ≤0.15 Å and non-Ag framework RMS ≤0.15 Å after a common framework translation is removed. Candidate roles were checked pairwise using the same rules. Contact checks compare each species-pair minimum with the smallest same-composition minimum in existing split/input geometries, allowing at most 0.05 Å compression; this is an empirical geometry guard, not a universal bond criterion.

## Role and lineage limits

All nine candidates are regenerated from original parent inputs. Each interface has one development-validation and one withheld-test geometry with interface-specific lateral Ag translations plus a +0.15 Å outward displacement of the marked framework atom. The second withheld-test candidate is a distance-axis shift of -0.20 Å (AgC), -0.175 Å (AgSi) or +0.175 Å (AgTi), measured along the documented marked-pair normal. Earlier wide-span candidate files provide the reference normal and lineage only; the output geometry is regenerated from the parent input. Per-ID displacement vectors, parent input hashes, original parent frame hashes, geometry hashes, contact minima and nearest inventory hits are in role_manifest.json.

These candidates are isolated from committed train/validation/test geometries and DFT input geometries under the stated thresholds, but they all inherit existing small-cluster Ti3SiC2 parent motifs. They do not establish morphology independence or thermal coverage. The six combined transverse/framework geometries are unrelaxed controlled designs. Geometry screening cannot establish DFT convergence or physical stability. An earlier V16 proposal set was excluded after each of its three files matched the corresponding frozen V16 blind input exactly; the geometry-only evidence is retained in attempt_diagnostics/attempt_01.

No DFT output, unseen reference label, MACE inference/training, MD or TTM was read or run. DFT labels require a separate owner assignment. Keep withheld-test labels sealed until model/checkpoint selection is frozen.

## Reproduction

Run the geometry-only generator with ASE and SciPy available:

    cd /workspace/-
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /workspace/.venvs/gpaw-mpi/bin/python -B research/mace-v12-transfer/coordination/reports/window-b/v17_split_geometry_design/final_set/generate_split_candidates.py --repo-root /workspace/- --output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v17_split_geometry_design/reproduction_01

The generator refuses an existing output directory, verifies source/manifest hashes and rechecks all nine geometry files, IDs, cell/PBC, duplicates and contact minima.
