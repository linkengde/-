# v13 force-error localization

Model SHA256: `85c16e3bf920418001377cf7049a788e201ac50811a4b2b639f40507b8bd5ae8`

Residuals are ranked by atom; the local region includes atoms within 4.5 Å of either marked contact atom.

| Label | All-atom vector RMSE (eV/Å) | Local RMSE (eV/Å) | Outside RMSE (eV/Å) | Max atom error (eV/Å) | Contact separation error (eV/Å) |
|---|---:|---:|---:|---:|---:|
| AgC_lateral_registry_holdout_v13 | 0.0891 | 0.1015 | 0.0604 | 0.2173 | 0.1749 |
| AgSi_lateral_registry_holdout_v13 | 0.0769 | 0.0795 | 0.0739 | 0.2734 | 0.1287 |
| AgTi_registry_probe_v13_01 | 0.0626 | 0.0712 | 0.0512 | 0.1281 | 0.0439 |
| AgC_registry_strain_acq_v14_01 | 0.1338 | 0.1606 | 0.0765 | 0.4167 | 0.3575 |
| AgSi_registry_strain_acq_v14_01 | 0.1021 | 0.1080 | 0.0944 | 0.2979 | 0.1618 |

Top-error atom IDs, elements, force directions and nearest neighbors are in the JSON. Use them to choose targeted follow-up DFT geometries; do not treat this narrow analysis as validation.
