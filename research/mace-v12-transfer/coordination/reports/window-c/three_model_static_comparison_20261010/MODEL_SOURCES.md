# Model source, license, and compatibility record

## Frozen inference environment

Existing analysis environment reused without reinstall: Python 3.12.14; PyTorch 2.5.1+cpu; mace-torch 0.3.16; ASE 3.29.0; NumPy 1.26.4; SciPy 1.17.1; CPU float64 inference. Both locally available models loaded and evaluated successfully. `r_max`, interaction count, species mapping and weight SHA256 are recorded in `model_inventory.csv`.

## MACE-MP-0b3-medium

- Official upstream checkpoint: <https://github.com/ACEsuit/mace-foundations/releases/download/mace_mp_0b3/mace-mp-0b3-medium.model>
- Project record: `research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v12_interface_energy/foundation_models/NOTICE.md` identifies the file as unmodified and reproduces the upstream release URL.
- Local model SHA256: `2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012`.
- The official MACE v0.3.16 catalog lists 89 elements, MPTrj training data, PBE+U, MACE >=0.3.10 compatibility and MIT license. Model is usable in the installed 0.3.16 runtime; its ordered `atomic_numbers` table contains C=6, Si=14, Ti=22, Ag=47.
- Loaded artifact reports 6.0 Å cutoff and two interaction layers. These values describe this exact local checkpoint, not every model in the MACE family.

## Frozen V18

- Project artifact: `research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v18_clean_core/training_clean_core_01/models/MACE_periodic_v18_clean_core.model`.
- Local SHA256: `757586671174db85acbf7971ccfe92de19cc58378ee2431a797efa6d6163cb62`.
- Archived as frozen epoch-79 project derivative. Loaded under MACE-Torch 0.3.16 on CPU float64; ordered species `[6, 14, 22, 47]`, 6.0 Å cutoff and two interaction layers.
- The project record did not establish an independent license for the derivative weights. This comparison does not infer that the MIT license on the MACE software or starting checkpoint automatically governs V18 weights.

## MACE-MATPES-PBE-0

- Official catalog source checked at MACE v0.3.16, upstream source commit `4d2da09413ac1407f37cdbb6b81fa28e4c15655e`: <https://github.com/ACEsuit/mace/blob/v0.3.16/README.md#latest-recommended-foundation-models>
- Catalog declares 89 elements, MATPES-PBE data, DFT(PBE) with no +U, a medium model, MACE >=0.3.10, and ASL licensing.
- Official checkpoint listed by that catalog: <https://github.com/ACEsuit/mace-foundations/releases/download/mace_matpes_0/MACE-matpes-pbe-omat-ft.model>
- The installed runtime satisfies the catalog's minimum version, but the checkpoint artifact was not retrieved: project/license scope has not been confirmed and the metadata endpoint attempted in the initial source check returned HTTP 403. No alternate host, mirror, or license bypass was used.
- Consequently, no artifact SHA, exact element order, cutoff, interaction count, or weight-level compatibility check exists, and no MATPES force/energy output is reported. Any future use requires the license scope to be resolved before download/inference and a new artifact hash/compatibility record.

The catalog describes the MACE software and linked model artifacts separately; software version compatibility does not settle model-weight licensing or artifact compatibility.
