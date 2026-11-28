# FVM_READINESS

## Stage decision

| Scope | Decision | Evidence and remaining condition |
|---|---|---|
|50³ QGS transient enthalpy regression|**GO for exploratory numerical work**|Same 120-step cycle and dt schedule; Numba/EW-safe agrees with reference liquid amount within 2.3×10⁻⁹ relative and normalized energy residual is below 10⁻⁹. Mesh convergence and physical inputs are separate.|
|100³ QGS 1 ns heat + 2 μs cool|**Conditional GO for performance development**|Completed in 553 s, 512 MiB RSS, 24 rejected trials, normalized energy residual 1.43×10⁻⁹. The 100³ earlier implementation exceeded 60 minutes without completion. Need time-step and grid convergence before paper-grade physical claims.|
|150³ QGS short stage|**GO only for bounded solver probes**|16 accepted steps in 172 s; 1.09 GiB RSS; reached 0.715 ns near melt onset; 5 rejected trials; checkpoint is restartable.|
|150³ full cycle|**NO-GO for this stage**|The short stage approaches the nonlinear melt band and dt drops to 2.34 ps. No end-to-end time/physical-progress evidence supports a full-cycle budget.|
|Ag/Ti₃SiC₂ real material prediction|**NO-GO**|The existing source/citation audit still lacks complete primary-source locations and Ag liquid data; preserve `MATERIAL_PROPERTY_AUDIT.md`.|
|Real arc heating / liquid flow / weld-resistance claim|**NO-GO**|No measured arc boundary, fluid solver, Marangoni or pressure terms, evaporation, or experiment coupling is included.|

## Why the GO is limited

The 100³ and 150³ results validate solver execution and energy accounting under candidate inputs. They do not establish that the selected local volumetric pulse represents an arc, that the 150³ geometry resolves every fine neck, or that transient differences are caused solely by skeleton size. 100³ EW vs EW-safe endpoint liquid fractions differ by about 0.56%, and there is no new mesh/time convergence for those two 100³ runs.

The time-step controller only retries on nonlinear failure and increases dt after convergence; a local truncation-error and liquid-fraction-change controller remains unimplemented. Likewise, an enthalpy-primary/active-set method has not yet been compared, and AMG was rejected as the default based on its small hot-state setup overhead. These are explicit follow-up items, not completed capabilities.

Next useful work is a controlled time-step study through melt onset, checkpoint I/O profiling, and a 150³ schedule only after melt-band cost is measured. Keep the solver's source mapping and material audit fixed during algorithm comparisons. Do not add molten-flow terms until the enthalpy solver is converged and material data is traceable.
