# Physical-model review: Ag–Ti₃SiC₂ single-temperature enthalpy study

Review date: 2026-10-28  
Repository baseline: `7035293a06a70c2515d3ab2c973417e48e5a1958`  
Scope: source definition and units, physical time-scale interpretation, candidate-property provenance, and the Case A/Case B limits. This is a code-and-record audit; no transient simulation was started.

## Main finding

The present calculation is a useful prescribed-source thermal-response model, but its source is not yet calibrated to a contact-closing event. It applies a fixed spatial heat-generation field directly to Ag enthalpy. It does not derive that field from a measured current/voltage waveform, contact resistance, arc energy, or an optical/electrical absorption model. In addition, the active Case A model sets Ti₃SiC₂ heat capacity and conductivity to zero. The current results therefore test heat localization in an Ag-only conducting network under an imposed load; they do not yet establish a quantitative melt threshold, a real contact cooling time, or a micro-versus-nano anti-welding ranking.

## Verified implementation facts

The solver in `phase_field_qsgs/qsgs_enthalpy_solver.py` defines

\[
q'''(x,y,z)=\phi_{Ag}(x,y,z)\frac{F_0}{D}
\exp\!\left[-\frac{(x-x_0)^2+(y-y_0)^2}{R^2}\right]
\exp\!\left[-\frac{z_{top}-z}{D}\right],
\]

where `F₀` is the input peak-flux-like amplitude, `R` is the radial Gaussian length, `D` is the depth e-folding length, and `φAg` is the Ag mask/fraction. For the native geometry, `qsgs_source_integration.py` analytically integrates this separable field over the binary Ag voxel mask. The solver then adds `q''' Δt` directly to the enthalpy-density field. Thus the source is a **volumetric generation term**, not a boundary-flux condition; the “flux” unit is used as its normalization before division by `D`.

The input is held constant during the heating interval and set to zero during source-off cooling. No absorptivity, reflectivity, current density, electric field, arc-plasma coupling, or source-partition parameter appears in this source definition. With the source mask, the reported “actual source power” is the integral of this prescribed field over Ag, not a measured device power.

The model has one temperature field. Heat is deposited directly into lattice enthalpy, and Ag uses sensible enthalpy plus a fusion-latent-heat plateau. It does not solve electron energy, electron–phonon transfer, melting-flow dynamics, vaporization, or decomposition. The native-grid outer boundaries are periodic in x/y and insulated in z.

## Source-unit and footprint audit

| Input / derived quantity | Value for the current study | Interpretation |
|---|---:|---|
| Peak amplitude `F₀` | 5 eV/(ps·Å²) = 8.01088317×10¹³ W/m² | Prescribed peak energy-flux scale; not an experimental measurement |
| Gaussian definition | `exp(−r²/R²)`, `R = 30 Å = 3 nm` | Radial e-fold length is 3 nm; half-maximum radius is 2.498 nm and FWHM diameter is 4.995 nm |
| Depth factor | `exp(−d/D)`, `D = 20 Å = 2 nm` | Depth e-fold length; the solver places the maximum at the top of the z-domain |
| Peak volumetric rate `F₀/D` | 4.005441585×10²² W/m³ | Before multiplication by the Gaussian, depth factor, and Ag mask |
| Peak fluence for a 1 ps flat pulse | 80.1088 J/m² | `F₀ × 1 ps`; a prescribed area-energy scale, not a calibrated contact dose |
| Infinite-plane, semi-infinite-depth integral | 2.265 mW | `F₀ πR²`; the ideal fully absorbing continuum value for this normalization |
| 70:30 geometry integrated source power | 1.551358946 mW | Geometry- and mask-integrated Ag source power recorded by the model |
| Energy at 0.2 ps / 1.0 ps for that geometry | 3.10272×10⁻¹⁶ J / 1.55136×10⁻¹⁵ J | Constant source power times duration |

For this Gaussian convention the FWHM diameter is `2 R √(ln 2)`, not `2R`. The exponential has no finite cutoff; the quoted `D` is an e-fold depth, not a hard absorption-layer thickness. The 70 vol% model integrates less power than the ideal 2.265 mW because the actual Ag mask and finite domain do not fill the whole source distribution. The mask is applied without renormalizing the field to a fixed total power.

### What is needed to connect this load to contact closing

For resistive contact heating, the eventual model input should be derived from an electrical/contact calculation, e.g. local `q''' = J·E` with temperature-dependent electrical properties and constriction geometry, or from a calibrated contact power history. The total resistive energy can be constrained by the measured contact voltage/current history, with the contact contribution separated from lead and arc contributions. Arc heating needs a separate spatial partition model or calibrated surface heat flux. Relevant experimental inputs include synchronized current and voltage, contact resistance/voltage drop, force and contact evolution, specimen geometry, arc occurrence/duration, and measured hot-spot dimensions.

The current 3 nm radius and 2 nm depth are therefore **numerical study inputs** until linked to such evidence. For a fair comparison, source footprint and time history must be the same in physical units in the micro/nano pair. A second, explicitly labelled equal-absorbed-energy comparison can normalize total Ag-deposited energy to isolate transport/storage effects; it must not be confused with the physical-load comparison.

## Single-temperature assumption and picosecond loading

The current 1 ps drive is applied straight to the Ag lattice. An electrical current initially transports energy through electronic degrees of freedom; for a sufficiently short load, electron and lattice temperatures need not be equal while energy is being deposited and redistributed. A single-temperature model is defensible only after showing that the relevant load and observation times are long compared with the material- and state-dependent electron–phonon equilibration time, and that local Fourier conduction is adequate at the source/feature length scale. The present repository does not supply or solve for electron heat capacity, electron–phonon coupling, electronic transport length, or an Ag-specific equilibration time, so it does not establish that criterion for a 1 ps pulse and a 3 nm source radius.

This is a validity question, not proof that the current result is wrong. If the actual contact heat pulse lasts much longer than a validated electron–lattice equilibration time, direct lattice heating may be an acceptable reduced model for the later response. If pulse and equilibration times are comparable, or if the source/gradient length is below the range where local Fourier transport is valid, compare against an electron/lattice two-temperature or nonlocal treatment. The existing independent double-temperature project is not modified by this review.

Allen’s electron–phonon relaxation treatment and the later nonequilibrium-metal model by Lin, Zhigilei, and Celli are relevant methodological references. They motivate material- and temperature-dependent relaxation checks; no numerical Ag relaxation time or coupling coefficient has been extracted from them in this audit, and they are not being used as sources for the candidate Ag property table.

## Candidate material properties: provenance remains open

The checked files are `phase_field_qsgs/thermophysical_candidates.csv`, `thermophysical_candidates.json`, and `thermophysical_model.py`. The CSV contains summary attributions (“NIST/CDC”, “NIST recommended data”, and “Barsoum et al.”), but none includes a full publication/edition, table or page, specimen condition, uncertainty, or archived source record. An outbound Crossref/DOI request from this environment returned HTTP 403, so the original records could not be independently opened. No exact property citation is invented here; the values below remain candidates.

| Material / input | Candidate currently in code | Implementation / unresolved source question |
|---|---:|---|
| Ag melting point | 1235.08 K | Thermodynamic anchor; exact standard/table and uncertainty not attached |
| Ag fusion enthalpy | 11,297 J/mol = 104,729.66 J/kg | Used as latent heat; original record and temperature convention missing |
| Ag solid density | 10,490 kg/m³ | Fixed for the full solid/liquid enthalpy-volume calculation; thermal expansion is omitted |
| Ag solid Cp | 25.4 J/(mol·K) at 300 K; 31.8 at Tm | Converted using 0.1078682 kg/mol and linearly interpolated between two anchors |
| Ag solid conductivity | 429 W/(m·K) at 300 K; 175 at Tm | Linearly interpolated; exact source and measurement convention missing |
| Ag liquid Cp / conductivity | 33.47 J/(mol·K) / 175 W/(m·K) at Tm | Single-anchor constants above Tm; liquid conductivity is explicitly tentative |
| Ag liquid density | 9,346 kg/m³ at Tm | Candidate is not used by the current solver; phase-volume change is not modeled |
| Ti₃SiC₂ density | 4,520 kg/m³ | Candidate does not distinguish theoretical/crystal density from porous specimen bulk density |
| Ti₃SiC₂ Cp | 110 J/(mol·K) at 300 K; 155 at 1473.15 K | Code conversion with 0.195708 kg/mol gives 562.06 and 792.00 J/(kg·K); interpolation is linear |
| Ti₃SiC₂ conductivity | 37 W/(m·K) at 300 K; 32 at 1473.15 K | Code is scalar/isotropic; source direction, specimen porosity, and uncertainty are not recorded |

Before paper-level use, attach the original source for each property with exact bibliographic details, location/table, units, temperature range, direction, specimen density/porosity, and uncertainty. Keep intrinsic density separate from composite bulk density. Add high-temperature data or a justified constitutive law instead of extrapolating the two-anchor curves. In particular, solid Ag conductivity is based on two anchors, liquid Cp/k on one anchor, and the current phase model has no vaporization or decomposition. The 70:30 QSGS volume fraction is not a mass ratio; with the current provisional densities it corresponds to about 84.4 wt% Ag.

## Case A and Case B interpretation

- **Case A** is an intentionally severe insulating-scaffold limit. The code sets Ti₃SiC₂ conductivity and sensible capacity to zero. The harmonic shared-face conductivity is then zero at every Ag/Ti₃SiC₂ face, so heat neither enters nor passes through the scaffold. It is not a realistic Ti₃SiC₂ property model.
- **Case B** includes the candidate Ti₃SiC₂ capacity and conductivity and assumes perfect thermal contact: there is one continuous temperature and no explicit interfacial temperature jump/resistance. It is the opposite interface limit, not a measured contact condition.
- Both current cases use periodic x/y and insulated z. Source-off evolution redistributes energy within an insulated domain; it does not represent heat loss to an electrode, fixture, surrounding material, radiation, or convection.
- Neither case includes a measured external heat sink, finite interface conductance, temperature-dependent electrical source, thermal expansion, stress, fluid flow, or weld-bond formation.

## Recommended model hierarchy before new long transients

1. Close the Ag/Ti₃SiC₂ source table and heat-source calibration described above.
2. Retain Case A only as a zero-scaffold lower-coupling bound and current Case B as a perfect-interface upper-coupling bound; explicitly keep their assumptions in each output.
3. Add a finite-interface model with separate Ag and Ti₃SiC₂ continuum temperatures and equal/opposite interfacial exchange `q'' = G_K (T_Ag − T_SC)`, where `G_K` is in W/(m²·K). For a cell-centered finite-volume face, use the series resistance
   \[
   G_{face}=\frac{A_{face}}{\delta_{Ag}/(2k_{Ag})+1/G_K+\delta_{SC}/(2k_{SC})},\qquad Q_{face}=G_{face}(T_{Ag}-T_{SC}).
   \]
   This permits an interface temperature jump while conserving equal-and-opposite energy. `G_K → ∞` recovers the perfect-contact face; a finite range must come from a source or measurement. Do not select `G_K` to force a desired melt result. This is a two-phase thermal model, not electron–phonon two-temperature coupling.
4. Add realistic outer heat paths and source-off sinks only when their material/contact geometry and boundary histories are specified. Check the combined Ag + Ti₃SiC₂ + boundary energy ledger.
5. Only then plan a short, checkpointed heat-window study. Keep source, mesh/time-step, and property sensitivity separate so each uncertainty is interpretable.

## References for the time-scale discussion

These references concern electron–phonon relaxation/nonequilibrium methodology; they are **not** property-source validation for the candidate table and no Ag numerical value is quoted from them here.

1. P. B. Allen, “Theory of thermal relaxation of electrons in metals,” *Physical Review Letters* 59, 1460–1463 (1987), DOI: [10.1103/PhysRevLett.59.1460](https://doi.org/10.1103/PhysRevLett.59.1460).
2. Z. Lin, L. V. Zhigilei, and V. Celli, “Electron-phonon coupling and electron heat capacity of metals under conditions of strong electron-phonon nonequilibrium,” *Physical Review B* 77, 075133 (2008), DOI: [10.1103/PhysRevB.77.075133](https://doi.org/10.1103/PhysRevB.77.075133).

