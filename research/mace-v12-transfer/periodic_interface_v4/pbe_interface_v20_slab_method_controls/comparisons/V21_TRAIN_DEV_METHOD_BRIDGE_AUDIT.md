# V21 train/development method-bridge audit

## Scope and evidence

This is a read-only audit of the approved V19 train36/dev3 method metadata, existing geometry/source audits, and the completed V21 fixed-geometry controls. It did not read training or development extxyz payloads, access test/sealed data, run GPAW or MACE, or alter dataset rows or labels.

The source-backed method audit reports PASS for 36 training rows and 3 development rows; its dataset-changed flag is `False`. The V21 controls use the exact same AgSi development input SHA256 `bca0a8fc1a2ef4e83eebe62d18f0b77480e33692be8886d9910dfd8622dd4b92` and are compared with `AgSi_v19_independent_dev_01`.

## Approved split method summary

| Split | Rows | XC / basis | Cutoff | k-points | Smearing | PBC |
|---|---:|---|---:|---:|---:|---|
| train36 | 36 | PBE / PW | 500 eV | Gamma (1×1×1) | 0.10 eV | [True, True, True] |
| dev3 | 3 | PBE / PW | 500 eV | 6×6×1 | 0.10 eV | [True, True, False] |

The shared settings are PBE, 500 eV plane-wave cutoff, and 0.10 eV Fermi smearing. The split differs in k-point sampling and boundary conditions: train36 uses Gamma and full three-dimensional periodicity; dev3 uses a 6×6×1 slab mesh and nonperiodic z with an xy dipole correction.

## V21 exact-development-geometry controls

| Control | SCF / comparison status | Native ΔE (meV/atom) | Free-energy Δ (meV/atom) | Force RMS Δ (eV/Å) | Marked Ag–Si projection Δ (eV/Å) |
|---|---|---:|---:|---:|---:|
| `AgSi_dev34_k6x6x1_sigma0p05_v21` | converged, 53 SCF iterations; `SENSITIVITY_EXCEEDS_BUDGET` | -0.04073 | -5.65683 | 0.01495 | 0.02619 |
| `AgSi_dev34_k6x6x1_sigma0p20_v21` | converged, 55 SCF iterations; `SENSITIVITY_EXCEEDS_BUDGET` | 0.60050 | 22.05047 | 0.03001 | -0.04723 |
| `AgSi_dev34_k4x4x1_sigma0p10_v21` | converged, 58 SCF iterations; `SENSITIVITY_EXCEEDS_BUDGET` | 0.64383 | 0.53755 | 0.01006 | -0.01574 |
| `AgSi_dev34_k8x8x1_sigma0p10_v21` | converged, 63 SCF iterations; `PASS_WITHIN_PROVISIONAL_BUDGET` | 0.18029 | 0.13687 | 0.00282 | -0.00174 |

The V21 controls vary slab smearing or in-plane k-mesh while keeping the approved AgSi development geometry fixed. They quantify sensitivity for this interface geometry only; they do not isolate a universal conversion between Gamma/full-periodic and slab/nonperiodic-z labels.

## Interpretation and next data boundary

The train/dev method distinction is coupled to different geometry classes and physical boundary conditions. It should not be treated as a same-geometry numerical offset or corrected with a scalar energy shift. The development labels are explicitly a slab calculation; the training rows are periodic bulk-like configurations. V21 controls inform the slab recipe but do not measure the effect of changing the same training geometry from Gamma/full periodicity to slab boundary conditions.

Keep the approved V19 train36/dev3 split and all existing labels unchanged. Keep every V20/V21 numerical control outside training and validation. Before acquiring new model data or training a successor, first complete the required source-verified Ti3SiC2 parent review independent of COD 9009647 and the B-side dataset/role audit; then register a separate data-acquisition task with an explicit method and split design. No further DFT or model training is authorized by this read-only audit.

## Source file SHA256

```text
2d9da5a13dee6bd3a1dd7d5940c1732616506184266f53d6497cd6e58cf4d6ec  research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v20_slab_method_controls/comparisons/approved_v19_smearing_consistency.json
9ce131494052cad0587bf4b70be6377403d664b93f92cc86c1c883c3e106eea3  research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_independent_dev/geometry_audit.json
7b52a8986fb58cf3faf840b81dfba3f9a81bb4898f10e25ef9185be468a5d2ca  research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v19_independent_dev/input_manifest.json
c5a0af1c1d4e877a21412aa0e3fc27f73ef676724751577e5946b0370b73a5e5  research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v20_slab_method_controls/comparisons/V21_EXACT_DEV_GEOMETRY_REVIEW.md
c404b0bc0af48955d710a19e314633e6d5ed2147f8eab74b67177edc5b828119  research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v20_slab_method_controls/comparisons/comparison_AgSi_dev34_k6x6x1_sigma0p05_v21__AgSi_v19_independent_dev_01.json
a298401341480d6198e6aec3da7f22d63e24048919006cd98a6d683ed4b49438  research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v20_slab_method_controls/comparisons/comparison_AgSi_dev34_k6x6x1_sigma0p20_v21__AgSi_v19_independent_dev_01.json
8f0b0472d326b0102b17ff08f21197f87717674549d6a9dd6b5ec9a922a404f9  research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v20_slab_method_controls/comparisons/comparison_AgSi_dev34_k4x4x1_sigma0p10_v21__AgSi_v19_independent_dev_01.json
eb2e00ed53159d38fedbe0367669ad36f239b960e45e6379829caf1e2851cdbb  research/mace-v12-transfer/periodic_interface_v4/pbe_interface_v20_slab_method_controls/comparisons/comparison_AgSi_dev34_k8x8x1_sigma0p10_v21__AgSi_v19_independent_dev_01.json
```
