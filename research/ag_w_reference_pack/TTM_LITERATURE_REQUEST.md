# Literature needed for exploratory Ag–W melt/quench TTM

The intended goal is a qualitative graduate-level simulation of laser-driven Ag–W contact melting, cooling, and resolidification. The Dai EFS potential is accepted as an exploratory model with known limitations; this checklist tracks what the user-supplied PDFs already cover and what remains. Do not fill missing entries with invented constants.

## Priority 1: tungsten electronic properties

**Sources and tables now supplied:** user-provided XTANT-3 tables for Ag and W `Ce(Te)`, `G(Te)`, and `ke(Te)`; Lin, Zhigilei & Celli (2008), including a candidate W curve digitized from Fig. 7; Daraszewicz et al. (2014) for an independent W low-`Te` coupling anchor and conductivity closure; and Medvedev & Milov (2020) as the source context for the XTANT coupling model. The input tables, upstream README/license, per-file hashes, and comparison are archived in `ttm_data/XTANT3_TTM_parameters_v3/` and `ttm_data/XTANT3_AG_W_DATA_AUDIT_20261124.md`.

The W `G(Te)` options materially disagree: at sampled common temperatures the supplied XTANT table is 3.8–17.9 times lower than the Lin digitization. The supplied standalone `Ce` values also differ from the separate `Ce` columns bundled with the K tables because they use different calculation setups. Keep these curves separate and compare explicitly; they still need a model choice and sensitivity run. The source K tables have only `Te` dependence, not `ke(Te,Tl)`.

For the remaining gaps, find temperature-dependent data or a documented model for:

- A `ke(Te,Tl)` relation that includes lattice-temperature dependence, or a justified temperature-only approximation for the exploratory scenario. The XTANT K tables and Duffy closure are now present but are not a validated high-nonequilibrium surface.
- Independent Ag `Ce/G` numeric data for cross-checking the supplied XTANT tables; Lin et al. plot Ag curves and cite tabulated data, but the linked university data host was inaccessible in this environment.
- An independent W `Ce(Te)` calculation or tabulated DOS result, if a stronger cross-check than the current plot digitization is desired.
- The temperature range, units, assumptions, and whether the values are measured, calculated, or fitted.

Search phrases for any additional W parameter work:

- `tungsten electron heat capacity electron thermal conductivity electron phonon coupling two temperature model`
- `laser heated tungsten electron phonon coupling G(T_e) thermal conductivity`
- `tungsten two-temperature model Ce ke G electronic temperature`

## Priority 2: Ag/W electron transfer across the interface

**No direct Ag/W interface conductance was identified in the supplied PDFs.** Li et al. (2024) adds measurements for Au/Pt and Au/Al and demonstrates strong interface/oxide sensitivity, but provides no Ag/W value. Look for Ag/W electronic thermal boundary conductance, electronic interface resistance, or a validated transport model that supplies an interfacial electron-transfer term. Record interface orientation, preparation/roughness, temperature range, and whether the reported conductance separates electronic and phononic channels. Giri et al. (2023) is useful background but is not a direct Ag/W value.

Search phrases:

- `Ag W interface electronic thermal conductance`
- `silver tungsten electronic thermal boundary conductance interface resistance`
- `Ag/W interface electron transport first principles thermal transport`

If direct Ag/W data cannot be found, provide papers that constrain a defensible range or a model for estimating it. The eventual calculation should compare a sensitivity range. Do not treat an arbitrary interface conductance as measured data.

## Priority 3: Ag reference data and laser pulse conditions

For Ag, retain the original citations and temperature ranges for `Ce`, `ke`, and `G` from the supplied XTANT tables and the separate Ag–Ti₃SiC₂ project, then check that the same definitions and units are used for W. The supplied XTANT tables are archived with a source notice that limits the dataset to non-commercial peaceful use.

For an experiment-specific comparison, find the laser wavelength, pulse duration and temporal shape, incident fluence/energy density, spot size, and optical reflectance/absorptance for Ag/W or a closely related contact. The current inherited spatial profile does not define `S(t)` or total pulse energy. If no matching experimental source exists, document the pulse as a chosen exploratory loading condition.

Search phrases:

- `Ag W contact laser pulse melting fluence`
- `silver tungsten laser irradiation pulse duration reflectivity`
- `Ag-W electrical contact laser melting resolidification`

## Priority 4: melt/quench checks for the EFS model

The pack contains W thermophysical data from Tolias et al. and the uploaded NASA TM 87568 program-task report (including its tungsten high-temperature measurement section), Ag thermochemistry from Pankratz, Jiao et al. on dilute nanocrystalline Ag–W, the calculated 1988 Ag–W phase diagram, Torre et al. (2021) on W–Ag segregation/demixing, and Waseda & Ohtani's liquid Ag structure data. The 1988 phase diagram explicitly rests on model assumptions and notes absent experimental thermodynamic data; it is not an experimental phase diagram. Additional liquid Ag–W mixing data would help interpret the model. The Dai EFS W melting point is about 21.8% above the reference; report any model melt temperature as model-specific unless calibrated.

Search phrases:

- `Ag W binary phase diagram liquid immiscibility mixing enthalpy`
- `silver tungsten phase diagram solubility`
- `tungsten melting latent heat liquid structure silver melting latent heat`

## What to send back

PDFs, supplementary files, DOI, or exact paper title are all useful. Prefer papers that expose tables, fitted functions, or enough methods to reproduce the values. For each paper, note which parameter it supports and the table/figure/page. We will record extracted values, unit conversions, uncertainty/range, and source hashes in the project before using them.

## Implementation limitation already identified

The stock LAMMPS `fix ttm/mod` in this environment does not support the requested `p p f` free-surface boundary. The existing finite-volume TTM code is Ag-only. Thermal-property tables alone do not complete the implementation; the solver must also support the free z surface, species-dependent parameters, the selected EFS model, and the chosen Ag/W interface closure, followed by an energy-ledger check. The full-size QSGS structure still fails static acceptance, so no temperature dynamics or laser run has started.
