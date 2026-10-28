# Contact-closing heat-source data audit

Audit date: 2026-10-28  
Scope: repository records for the current Ag 70:30 volume-fraction thermal model. No new source history or transient was generated.

## Data found in the current model

The thermal solver’s source is a prescribed, constant-during-heating field with the spatial form

\[
q'''=\phi_{Ag}\frac{F_0}{D}\exp(-r^2/R^2)\exp(-d/D),
\]

using `F₀=5 eV/(ps·Å²) = 8.01088317×10¹³ W/m²`, `R=30 Å=3 nm`, and `D=20 Å=2 nm`. It is switched off at the requested heat duration. For the saved Ag 69.9996 vol% mask, the model integrates this to 1.551358946 mW. This is a numerical heat source delivered directly to Ag enthalpy. The “absorbed power” label describes the integral of that imposed field over the Ag mask; it is not experimentally measured contact power and does not model reflection/absorption or current flow.

The 70:30 QSGS run files record heat duration, mask, and power. In the checked model/report directories, no synchronized closing-current waveform, contact voltage drop, contact-resistance trace, arc voltage/energy, or thermal calibration of the 3 nm × 2 nm source footprint was found. The code parameter set does not contain these measurements.

## Required experimental/electrical inputs

| Input | Why needed | Status in checked model/report tree |
|---|---|---|
| Synchronized `I(t)` and voltage directly across the contact during closing | Constrain total event energy and separate contact voltage from source/leads | Not found |
| Contact resistance or contact voltage `V_c(t)` and its measurement bandwidth | Estimate `P_contact(t)=I(t)V_c(t)`; resolve transient constriction/contact changes | Not found |
| Arc voltage/current and arc duration, if an arc occurs | Separate arc deposition from constriction Joule heating | Not found |
| Switch/contact force, closing velocity/displacement, electrode geometry and material | Constrain contact area, constriction geometry, and load history | Not found |
| Contact surface roughness, oxide/coating state, and contact spot/hot-zone size | Constrain local current concentration and heat partition | Not found |
| Temperature-dependent electrical resistivity/conductivity of Ag and Ti₃SiC₂ | Compute `q'''=J·E` and its phase allocation | Not source-verified in this model tree |
| Energy partition to the two contacts, bulk electrodes, arc/plasma, and fixture | Map measured electrical work to the local microstructure domain | Not found |
| Independent temperature/spot constraint (high-speed imaging, calibrated pyrometry, or other suitable diagnostic) | Check inferred spatial heat-source footprint and peak response | Not found |

## Recommended source construction

1. Acquire synchronized current and voltage traces, with separate contact voltage measurement and a parasitic/lead calibration. Integrate `E_contact=∫ I(t)V_c(t)dt`; report the measurement bandwidth, offsets, and uncertainty.
2. Record mechanical contact timing/force and electrode dimensions. Use these to create the time-dependent constriction geometry and electrical boundary condition.
3. Solve current continuity with phase- and temperature-dependent electrical properties, then map local Joule heating as `q'''=J·E` to the resolved Ag and Ti₃SiC₂ geometry. Keep external contact and bulk-electrode losses in a coupled macro model or explicit boundary condition.
4. If an arc contributes, measure its voltage/current/time separately and use a documented surface heat-partition model. Do not combine arc energy with contact Joule energy under one arbitrary Gaussian.
5. Calibrate spatial deposition against a measured hot-spot/contact scale and energy balance. Only after this should a reduced Gaussian/exponential source be fitted; preserve its integrated energy and uncertainty.

## Interim rule

Keep `F₀`, `R`, `D`, and the 1 ps pulse labeled exploratory. Do not describe them as a real contact-closing load, infer an experimental melt threshold from them, or adjust them to force a melt. Until the required traces exist, a sensitivity matrix over clearly labeled source inputs is useful for software verification but cannot establish a physical weld ranking.

