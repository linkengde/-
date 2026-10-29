# Ag70 1× Case B model preflight

## Outcome

**CASE_B_PREFLIGHT_PASSED** on the native Ag70 QSGS geometry. This validates the existing Case B implementation for short, source-off conduction tests; it does not establish a melt-pool result.

## Model used

- Geometry: `runs/single_temperature_mesh2_validation_20261027/qsgs_morphology_ag70vol.npz`; SHA-256 `03ef47b1dd534b3215cf03f42c8cf7e887dcdfcc185b33bab7453cd838c1faac`.
- Grid: `[100, 100, 50]`; mesh refinement `1×`.
- Ag/Ti₃SiC₂: `69.999600%` / `30.000400%` by actual grid count.
- Initial state: 300 K; x/y periodic; z insulated.
- Case B: both phases store and conduct heat with ideal thermal contact; only Ag receives the source.
- Source: peak flux `8.01088e+13 W/m²`, Gaussian radius `30.0 Å`, absorption depth `20.0 Å`; integrated power `0.00155135895 W`.
- Explicit diffusion limit at 300 K: `0.0404211 fs`.

## Checks

| Check | Result |
| --- | --- |
| Native shape, Ag70 geometry identity, 1× resolution | passed |
| Heat source has no support in Ti₃SiC₂ cells | passed |
| 300 K, no source, 0.001 ps static run | passed; energy `0.0000e+00 J` |
| Uniform 500 K, no source, 0.001 ps | passed; max drift `5.684e-14 K` |
| Hot Ag / cold scaffold transport | passed; scaffold max `300.000000 → 332.040304 K` |
| Hot Ag cools while scaffold stores heat | passed; Ag mean `600.000000 → 595.653614 K`; scaffold gain `1.1632e-16 J` |
| Closed-system energy conservation in conduction test | passed; relative change `2.038e-16` |
| Explicit timestep guard | passed |
| Ti₃SiC₂ upper-property limit | passed; states above 1473.15 K rejected |

## Material and physical limits

The solver uses the repository's first-pass Ag and Ti₃SiC₂ candidate properties. The Ti₃SiC₂ table ends at 1473.15 K, and the model rejects states above that limit. These candidates are not fully source-verified. Case B assumes perfect Ag/Ti₃SiC₂ contact and an insulated outer boundary; it is not a simulation of heat loss to the electrodes or ambient. A shared single temperature, no liquid flow, no vaporization, no expansion, and no mechanical constraint are included.

The next numerical gate is a 0.01–0.05 ps Ag-only source-on run from a fresh 300 K Case B state. Longer segmented heating and a source-off cooling branch should proceed only while the temperature/property limits and energy ledger remain valid.
