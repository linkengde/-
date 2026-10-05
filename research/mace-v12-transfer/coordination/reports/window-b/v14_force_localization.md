# v14 force-error localization

Model SHA256: `41bd411355b034b0dab85dbdbda932ecfc394122cb67a2b51a550bd04207645d`

Post-hoc diagnosis of the published blind-holdout screen; local means within 4.5 Å of either marked contact atom.

| Label | All-atom vector RMSE (eV/Å) | Local RMSE (eV/Å) | Outside RMSE (eV/Å) | Max atom error (eV/Å) | Contact separation error (eV/Å) |
|---|---:|---:|---:|---:|---:|
| AgC_registry_holdout_v14_01 | 0.1537 | 0.1762 | 0.0937 | 0.3241 | 0.0614 |
| AgSi_registry_holdout_v14_01 | 0.1178 | 0.1266 | 0.1079 | 0.2793 | 0.2349 |
| AgTi_registry_holdout_v14_01 | 0.1466 | 0.1775 | 0.0834 | 0.4652 | 0.1926 |

Top-error atom IDs, force directions and neighbors are in the JSON. Use this only to guide follow-up DFT; the v14 screen remains failed.
