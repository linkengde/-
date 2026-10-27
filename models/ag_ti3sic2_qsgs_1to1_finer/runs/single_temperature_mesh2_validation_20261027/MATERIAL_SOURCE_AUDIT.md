# Thermophysical source audit for the Ag 70:30 vol% screening runs

Audit date: 2026-10-27 UTC. Values below are the candidates currently consumed by `phase_field_qsgs/thermophysical_model.py`; this audit does **not** certify them as literature-validated. The repository has summary attributions but no exact publication/table/page or archived source record. A direct Crossref API lookup was attempted from this execution environment and the network proxy returned HTTP 403, so source records could not be fetched or checked here. No bibliographic details have been guessed.

| Material / quantity | Candidate value(s) and units | Stated attribution in current files | Status / model treatment |
|---|---:|---|---|
| Ag molar mass | 0.1078682 kg/mol | Candidate config | Verify against a named standard; used for molar-to-mass Cp conversion. |
| Ag melting point | 1235.08 K | User-provided summary attributes thermodynamic anchors to a NIST/CDC table | Exact table/version missing; candidate, not independently verified. |
| Ag fusion enthalpy | 11,297 J/mol; 104,729.66 J/kg | Same general NIST/CDC attribution | Exact record missing; used as the isothermal enthalpy plateau. |
| Ag reference density | 10,490 kg/m³ | Same general NIST/CDC attribution | Temperature/source convention not archived; solver holds this solid density fixed. |
| Ag liquid density at melting | 9,346 kg/m³ | User-supplied candidate | Provenance missing; not used in the thermal solver's volume capacity, but relevant to conditional mass/volume conversion. |
| Ag solid Cp | 25.4 J/(mol K) at 300 K; 31.8 at 1235.08 K | Same general NIST/CDC attribution | Only two anchors; linearly interpolated from 300 K to Tm. Exact records missing. |
| Ag solid conductivity | 429 W/(m K) at 300 K; 175 at 1235.08 K | User summary says NIST recommended data | Exact table/temperature convention missing; linearly interpolated. |
| Ag liquid Cp | 33.47 J/(mol K) at Tm | User-supplied candidate | Single anchor; solver holds it constant above Tm. Temperature dependence is unknown. |
| Ag liquid conductivity | 175 W/(m K) at Tm | User summary calls this a tentative recommendation | Single tentative anchor, held constant above Tm; high uncertainty and no verified citation. |
| Ti₃SiC₂ molar mass | 0.195708 kg/mol | Formula-derived candidate | Reproducible from formula atomic weights; document atomic-weight standard in final data package. |
| Ti₃SiC₂ density | 4,520 kg/m³ | User-supplied candidate | Exact specimen density/porosity and source not archived. |
| Ti₃SiC₂ Cp | 110 J/(mol K) at 300 K; 155 at 1473.15 K | User-provided summary attributes values to “Barsoum et al.” | Citation, specimen, phase purity and table/page missing. Linear interpolation is a model assumption. Dividing by 195.708 g/mol gives 562.06 and 792.00 J/(kg K); the earlier 656–924 conversion is inconsistent. |
| Ti₃SiC₂ conductivity | 37 W/(m K) at 300 K; 32 at 1473.15 K | User-provided summary attributes values to “Barsoum et al.” | Exact citation and directional/porosity convention missing. Linear interpolation is a model assumption. |

## Applied ranges and consequences

- The Ag solid Cp/k anchors cover only 300 K to the configured melting point. The liquid Cp/k have only one anchor at melting and are treated as constants. The enthalpy model omits vaporization and decomposition.
- Ti₃SiC₂ anchors cover 300–1473.15 K. The existing code rejects extrapolation above that interval, but the requested Case A runs set the scaffold to zero heat capacity and zero conductivity, so those runs do not validate scaffold thermal response.
- The 70:30 volume mask converts to 5.415:1 Ag:Ti₃SiC₂ mass under the current 10,490/4,520 kg/m³ candidate densities (84.412 wt% Ag). A strict 1:1 mass composition would instead be about 30.113 vol% Ag under those same provisional densities.
- Heat-source values (5 eV/(ps·Å²), 30 Å Gaussian radius, 20 Å absorption depth) are study inputs, not measured contact-joule calibration. They remain unchanged for this study.

## Required closure before quantitative paper claims

Attach the exact primary source records (DOI/edition, table/page, unit and temperature), quote data uncertainties, verify solid/liquid density convention and porosity, and replace the generic “NIST/CDC” and “Barsoum et al.” labels with precise citations. Then re-run a sensitivity study for the parameters that control onset and liquid lifetime, especially source energy, Ag latent heat/Cp/k and the scaffold/interface thermal path. Until then all temperatures, phase onset and timings remain exploratory model outputs.
