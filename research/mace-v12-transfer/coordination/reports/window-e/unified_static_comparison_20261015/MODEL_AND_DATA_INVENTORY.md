# Frozen model and reference-data inventory

Run date: 2026-10-15 UTC. Branch: `cloud-e-b1-dft`. This report uses frozen weights only. No DFT, model training, potential modification, MD, or TTM was performed.

## Retrieved archives and structure verification

| Archive | Git commit and repository path | Archive SHA-256 | Verification |
| --- | --- | --- | --- |
| 77-atom three-gap package | `3c8ffa28cf6f5aaa72981d8b9f1a29389f9d8d50`; `research/mace-v12-transfer/coordination/reports/window-c/minimal77_ag_tisc_dft_input_preparation_20261014/MACE_SMALL_INTERFACE_COMPARISON_DATA_20261015.zip` | `f37d5c84b919ac583c9811b557c8eacdb06254d2d4095b6634c0f9ab927ba678` | SHA-256 and internal `SHA256SUMS.txt` pass; unpacked separately under `/tmp` |
| 135-atom DFT handoff | `1353d3efa04c1bdd19f139b3535f21d02985013d`; `research/mace-v12-transfer/coordination/reports/window-e/AG_TISC_C_DFT_HANDOFF_20261015.tar.gz` | `4e35ce85e6a47295e07c6e4e688925bc343fee7cdd966d39e3619b9d00d7026f` | SHA-256 and all 208 package entries in `SHA256SUMS.txt` pass; unpacked separately under `/tmp` |

The 77-atom files have 77 atoms each, composition Ag₂₈Ti₂₁Si₁₄C₁₄, identical cells/PBC/order, and unique `atom_id`s. The cell is `[[7.44683267,-2.57965851,0],[-1.48936653,7.73897552,0],[0,0,37.9451596]] Å`, PBC `TTT`. Their SHAs are 2.16 Å `51223cc8f1dd881a6fa1f2e4ea5564f886ca3682e6d4ef4389393d3a6768c5d3`, 2.40 Å `1de937b99b2d1f5b4be3b9b45e4d13e4e634bd99531110a080abf7ce970319dd`, and 2.64 Å `52815dc6d9f4868cbdd988f11fd86cb2fac47b7dc7faf56b2449887f88a93a34`. Only the Ag z coordinates move, by +0.24 Å and +0.48 Å relative to 2.16 Å; Ti₃SiC₂ coordinates and cell are unchanged. This is a single 0° strained parent family (TSC in-plane about −3.0045%, Ag about +3.0975%), not the 12.51983° low-mismatch 135-atom model. No 77-atom DFT labels were present.

The 135-atom relaxed geometry is Ag₂₁Ti₅₇Si₁₉C₃₈, PBC `TTT`, cell `[[13.38617866,0,0],[6.69308933,11.59277078,0],[0,0,37.9451596]] Å`. Optimized structure SHA-256: `eda2586c63c6c3c793ffce99417e7c04569941b0a2a916fda1b5402b7ba52d9f`; force table SHA: `efd00ea38c50e6bf9eb54ac8f45d7f0ed9bd8a2a31a8f31d725af7f74b7f03c6`. Full DFT converged at the final geometry (22 SCF steps, Fmax 0.0427799 eV/Å). Ag-only (21 atoms) and TSC-only (114 atoms) fixed-geometry DFT references also report converged; TSC-only was recovered from its completed checkpoint/log and has 44 SCF steps. All three use GPAW 26.7.0, ASE 3.29.0, PBE/PAW, 500 eV, Gamma 1×1×1, Fermi–Dirac 0.10 eV, spin-unpolarized, fixed cell and coordinates, and 1e-5 energy/density/eigenstate thresholds. Component coordinates/cell/PBC were checked against the full geometry. Source atom ID is the 1-based full-structure row ID; the extxyz `atom_id` field is retained separately in the per-atom CSV.

DFT total-energy and force-consistent free-energy references are kept distinct. Composition-subtracted increments are `Eint,total = -22.229525474612 eV` and `Eint,free = -21.378661277249 eV`. The model energy is compared to the native total-energy increment; the free-energy increment is separately shown as a convention sensitivity, not silently mixed. PAW setup hashes from the DFT run manifest: Ag `bff00ac79a7fa6e8c3a71c8b34d3a9c488e6bf1c00b3e9c4c91c231551145e97`, Ti `bf822f91317bebbb6ea9b456dc4e8f72fb2a8b824984221578bf1bd7048add2f`, Si `2defcf97f4072f0946dd4ad44cbb78bb48987a64f38a79afbf62ff66c6154ef0`, C `fbef4c5328af74d9fb1b5b7f9aa024cd9ec5a2b9b862351f5cff280e0f83c983`.

## Frozen model weights and data roles

| Candidate | Inference weight SHA-256 | Existing role/evaluation evidence |
| --- | --- | --- |
| MACE-MP-0b3-medium | `2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012` | Foundation model; no project-specific train/test split claimed |
| MACE V12 | `011280f8d999a62c1603848226fd64c9592740d54e28d4d04f76e2035a034628` | 23 train / 2 valid / 3 test; existing overall validation remains `UNDETERMINED`; small frozen-test results are not broad interface validation |
| MACE V14 | `41bd411355b034b0dab85dbdbda932ecfc394122cb67a2b51a550bd04207645d` | Existing assessment’s evaluated `run-45.model` (used here to match its recorded evaluation); assessment also names a different selected epoch-71 checkpoint SHA `d639c13499d86418d66f59b9be4800a3758bfdb68aeccdf9c60a0ecca27a9fa9`. Prior assessment is `SCREEN_FAIL` on interface force/separating-force gates |
| MACE V17 | `559b12c024fe8e2155e26ffc7ca48cf9fa1c2682ea37ff9c5367af14a5cc43d3` | 43 train / 3 reused development frames / no test file; provisional screening and source/method caveats remain; previous development gates failed |
| MACE V18 | `757586671174db85acbf7971ccfe92de19cc58378ee2431a797efa6d6163cb62` | 30 train / 3 development / no test; previous development gates failed; training set had no pure Ti₃SiC₂ frames |

Inference used ASE 3.29.0, MACE 0.3.16, Torch 2.5.1+cpu, CPU, float64, one Torch thread, ASE eV and Å for all candidates. All five weights loaded and produced finite energies/forces on all evaluated structures. Exact-geometry scans of V12/V14/V17/V18 `train`, `valid`, and available `test` extxyz partitions found no exact matches for the 135 full structure, Ag/TSC subsets, or the three 77 structures. This is an exact-geometry check only: a formal parent-family join to every historical scored development group is not present in all model manifests, so family-level independence remains **UNVERIFIED**. Do not call this a blind test.

The retained V18 pure-backbone warning is material: a separate 48-atom Ti₃SiC₂ DFT reference had MACE-MP force RMSE 0.038173 eV/Å and V18 2.888854 eV/Å (about 75.7×), with no pure TSC frames in V18's 30-frame training set. That independent existing result remains in force; the local interface response cannot override it.

## Static inference outputs

The 135-atom model table and per-atom forces are in the two CSVs; the 77-atom response table is the third CSV. Force vector RMSE is `sqrt(mean_i(||Fmodel_i−FDFT_i||²))`, in eV/Å; element errors use only that species. Phase and effective interaction forces are mapped by full-source atom ID to the same positions and cell. The DFT effective cross force is `Ffull−FAg-only` on Ag and `Ffull−FTSC-only` on Ti/Si/C; MACE uses the analogous component subtraction. This is an electronic-response increment, not a strict two-body decomposition.
