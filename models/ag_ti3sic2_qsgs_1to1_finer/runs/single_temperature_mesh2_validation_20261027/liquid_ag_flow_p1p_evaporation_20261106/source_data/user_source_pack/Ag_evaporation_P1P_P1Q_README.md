# Ag evaporation / recoil-pressure source pack for P1P-P1Q

## Included structured data
- ag_vapor_pressure.csv: Wang et al. (2024) liquid-Ag fit plus NIST Antoine comparison over overlapping ranges.
- ag_vaporization_enthalpy_derived.csv: derived from the Wang et al. vapor-pressure fit using Clausius-Clapeyron and ideal-vapor approximation; not direct calorimetry.
- ag_evaporation_coefficient_table8.csv: alpha_ev values reported in Table 8 (1373-1773 K).
- ag_measured_evaporation_rate_table5.csv: Table 5 measured rates at 1, 10, 50, 100 Pa.
- ag_evaporation_logistic_fit_table7.csv: Table 7 logistic fit coefficients.
- ag_vacuum_evaporation_kinetics.csv: overall mass-transfer coefficients and activation energy.

## Primary vapor-pressure correlations
Wang et al. (2024), liquid Ag, 1234-2436 K:
log10(p_sat/Pa) = -14400/T - 0.85 log10(T) + 13.825

NIST Antoine fit, 1823-2425 K:
log10(P/bar) = 1.95303 - 2505.533/(T - 1194.947)
NIST source uses C=-1194.947 K.

## Important interpretation
- Wang et al. (2024) is a pure-Ag vacuum-distillation study. Its actual rate and alpha values are not automatically transferable to an electrical-arc environment.
- alpha_ev values are 0.2118, 0.0834, 0.1428, 0.1329, 0.2715 at 1373, 1473, 1573, 1673, 1773 K. Do not extrapolate above 1773 K as validated data.
- k_Ag is an overall vacuum mass-transfer coefficient (m/s), not the Hertz-Knudsen alpha.
- Ag vapor partial pressure must be specified; do not substitute total shielding-gas pressure.
- Table 5 rate units are transcribed as g cm^-2 min^-1 from the article context; verify the publisher PDF table header before code ingestion.
- NIST page gives 1823-2425 K Antoine range and cites Hartmann & Schneider (1929). Inspect the underlying primary reference before treating this as independent confirmation.
- Derived vaporization enthalpy assumes ideal vapor and uses the slope of the Wang et al. vapor-pressure fit; label it derived, not directly measured.

## P1P/P1Q numerical cautions
Much et al. (2024) shows that diffuse-interface/CSF interface temperatures and temperature-dependent fluxes can be significantly biased by property jumps; it recommends careful interface-midplane temperature evaluation and improved parameter-scaled CSF treatment. This is particularly important because Ag vapor pressure depends strongly on surface temperature.
Wang et al. (2020) is a general metal evaporation/recoil-pressure model, not Ag-specific; its abstract reports limitations of conventional Anisimov treatment in near-vacuum.
Holyst & Litniewski (2009) is an MD study of liquid-film evaporation into vacuum; the PubMed abstract is available, but a full-text PDF was not verified as openly accessible.

## Official source/download links
1. Wang et al. (2024), pure Ag evaporation (open-access publisher landing page; choose View PDF):
https://doi.org/10.1016/j.nxmate.2024.100189
2. Much et al. (2024), open-access PDF:
https://link.springer.com/content/pdf/10.1186/s40323-024-00270-6.pdf
3. Wang, Zhang & Yan (2020), APS PDF endpoint:
https://journals.aps.org/prapplied/pdf/10.1103/PhysRevApplied.14.064039
4. NIST Silver Antoine data:
https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440224&Mask=4&Plot=on&Type=ANTOINE
5. Holyst & Litniewski (2009), PubMed record:
https://pubmed.ncbi.nlm.nih.gov/19239309/

## Download limitation
The ZIP contains structured numerical tables, equations, source metadata and official publisher links. It does NOT contain publisher PDFs or separate Supporting Information files. Those full external PDFs were not copied into this ZIP.
