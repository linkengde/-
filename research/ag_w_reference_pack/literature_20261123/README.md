# Ag–W laser/TTM/thermal transport literature intake

Research date: 2026-11-23. This is a source and access register, not a claim that all linked PDF binaries have been mirrored into this repository. Free full-text PDF links were opened and bibliographic metadata checked. The current execution environment could not retrieve external PDF binary bytes for upload, so the open PDFs are linked below; no false "downloaded to GitHub" status is assigned. Existing user-supplied Dai 2006, He 2018 and Jiao–Schuh 2018 PDFs already archived in `papers/` remain unchanged.

## 1. Tungsten two-temperature model (highest priority)

### Lin, Zhigilei & Celli (2008)
**Electron-phonon coupling and electron heat capacity of metals under conditions of strong electron-phonon nonequilibrium.** Physical Review B 77, 075133 (2008).
- Free full text PDF (APS Harvest): https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.77.075133/fulltext
- DOI: https://doi.org/10.1103/PhysRevB.77.075133
- Covers W and Ag; (C_e(T_e)), (G(T_e)), and a Drude (k_e(T_e,T_l)) form. W (C_e) begins to depart from linearity above ~3000 K, peaks in its deviation near ~11,000 K, and begins to saturate/decrease around ~17,000 K. The paper gives (gamma_{W}=137.3) J m⁻³ K⁻² for the low-(T_e) linear reference. It notes significant uncertainty/spread in W (G) values and a rise with (T_e).

### Daraszewicz et al. (2014)
**Determination of the electron-phonon coupling constant in tungsten.** Applied Physics Letters 105, 023112 (2014).
- Free author-hosted PDF: https://discovery.ucl.ac.uk/id/eprint/1451162/1/Duffy_1.4890413.pdf
- DOI: https://doi.org/10.1063/1.4890413
- W(100) pump-probe study. Reports (G_0=1.65\times10^{17}) W m⁻³ K⁻¹ from first principles and (1.4(3)\times10^{17}) W m⁻³ K⁻¹ from fitting reflectivity dynamics. Pump duration 90 fs; incident fluence used for the fit 5.22 mJ cm⁻². Its W conductivity model uses (k_e/k=0.75) at 300 K, total (k=174) W m⁻¹ K⁻¹, and (n=1.5) in the Drude-like (T_l) dependence.

### Medvedev & Milov (2020)
**Electron-phonon coupling in metals at high electronic temperatures.** Physical Review B 102, 064302 (2020).
- Free full text PDF: https://ris.utwente.nl/ws/files/222992584/Medvedev_2020_Electron_phonon_coupling_in_metals_.pdf
- arXiv: https://arxiv.org/abs/2005.05186
- DOI: https://doi.org/10.1103/PhysRevB.102.064302
- Includes W and Ag (G(T_e)) calculations using XTANT-3 dynamical coupling. It explicitly compares with Lin et al. and notes substantial model/experimental spread. Use as a distinct model branch, not a replacement for the 2014 W pump-probe calibration.

### Akhmetov, Milov & Medvedev (2024)
**Electronic heat conductivity in a two-temperature state.** International Journal of Heat and Mass Transfer 228, 125674 (2024).
- Free arXiv PDF: https://arxiv.org/pdf/2311.12458
- arXiv: https://arxiv.org/abs/2311.12458
- DOI: https://doi.org/10.1016/j.ijheatmasstransfer.2024.125674
- Calculates electronic conductivity to (T_e=40,000) K for W and Ag. The associated XTANT-3 repository (K_W.dat) and (K_Ag.dat) extend to 50,000 K. These are (T_e)-indexed data, not a full (k_e(T_e,T_l)) surface.

### XTANT-3 data release
- Public data repository: https://github.com/N-Medvedev/XTANT-3_coupling_data
- Zenodo record: https://zenodo.org/records/10183636
- DOI shown on Zenodo record: 10.5281/zenodo.10183636
- Contains separate Ag/W (C_e), (G), and (k_e) tables. Imported copies and source blob hashes are in `ttm_data/XTANT3_TTM_parameters_v3/`.
- Caveat: Zenodo explicitly warns room-temperature values may not always be reliable. Separate `Ce_*.txt` and the (C_e) column in `K_*.dat` use different calculations.

## 2. Ag/W interface electronic thermal transport

### Giri et al. (2023)
**Ultrafast and Nanoscale Energy Transduction Mechanisms and Coupled Thermal Transport across Interfaces.** ACS Nano 17, 14253–14282 (2023).
- Free author-hosted PDF: https://patrickehopkins.com/wp-content/uploads/2024/02/giri2023aa.pdf
- Open PMC record: https://pmc.ncbi.nlm.nih.gov/articles/PMC10416573/
- DOI: https://doi.org/10.1021/acsnano.3c02417
- Critical finding: the review says electronic metal/metal boundary conductance has not been studied for metals with relatively high phonon heat-conduction contributions, explicitly including W. Thus no direct Ag/W (h_{ee}) value was identified in this search.

### Electronic interfacial thermal conductance between metal nano-films (2024)
**Electronic interfacial thermal conductance between metal nano-films in sub-picosecond non-equilibrium transport process.** International Journal of Heat and Mass Transfer 229, 125717 (2024).
- DOI: https://doi.org/10.1016/j.ijheatmasstransfer.2024.125717
- No free full-text PDF was verified in this search.
- Analog measurements, not Ag/W: Au/Pt (h_{ee}=4.4\pm1.8) GW m⁻² K⁻¹; Au/Al (0.13\pm0.05) GW m⁻² K⁻¹, with the oxide layer strongly suppressing electronic transport.
- These values can motivate a sensitivity sweep only; they must not be labeled as Ag/W measurements. A provisional exploratory sweep could test 0.1, 1, 4.4 and 10 GW m⁻² K⁻¹; 10 is an extreme high-(T_e) theory-scale stress test, not a measured Ag/W value.

### Existing Ag/W interface DFT reference
- He et al. (2018), DOI: https://doi.org/10.1021/acs.jpcc.8b05920
- The user-supplied author PDF is already archived in `papers/He_2018_AgW_Interface_Stability.pdf`.
- It provides interface structure, work of separation and interface energy, but not an electronic thermal boundary conductance.

## 3. Laser conditions: closest free experimental analogues

### Byskov — ultrashort-pulse laser ablation of Cu, Ag and W
- Free PhD thesis PDF: https://phys.au.dk/fileadmin/site_files/publikationer/phd/jeppe_byskov.pdf
- Chapter 2 is based on the paper *Ultra-short pulse laser ablation of metals: threshold fluence, incubation coefficient and ablation rates*, Applied Physics A; DOI: https://doi.org/10.1007/s00339-010-5766-1
- Experiment: 100 fs, 800 nm Ti:sapphire pulses; 1/e spot size ~40 µm; pulse energies 9–185 µJ; repetition rate 10 Hz–1 kHz.
- Initial polished-surface Fresnel reflectance at 800 nm: Ag 0.981, W 0.496.
- Single-shot threshold fluences: Ag (1.5\pm0.4) J cm⁻²; W (0.44\pm0.02) J cm⁻².
- The fitted ~100-shot low-fluence thresholds are Ag (0.83\pm0.02) J cm⁻² and W (0.22\pm0.01) J cm⁻². Effective penetration depths from that ablation analysis are Ag (100\pm6) nm and W (25\pm1) nm; low-fluence effective absorptances are 0.294 and 0.960 respectively.
- Important: these are pure-metal ablation data with incubation/surface evolution, not laser welding of an Ag–W composite. The initial reflectance and multi-shot effective absorptance are different physical conditions.

### Daraszewicz et al. (2014) W pump-probe
- Free PDF and DOI above.
- 90 fs pump; 400 nm pump wavelength (12.5 nm optical penetration depth in W); fit uses 5.22 mJ cm⁻² incident fluence; reflectivity probes at 875 and 1000 nm.
- This is a W ultrafast relaxation experiment, not an Ag/W contact-welding experiment.

**No directly comparable Ag–W electrical-contact laser-welding experiment with verified pulse parameters was found in this search.** Treat the existing Ag-only spatial source as exploratory until a target experiment and pulse envelope are specified.

## 4. Ag–W phase stability, mixing enthalpy and melt references

### Vijayakumar, Sriramamurthy & Nagender Naidu (1988)
**Calculated phase diagrams of Cu–W, Ag–W and Au–W binary systems.** Calphad 12(2), 177–184 (1988).
- DOI: https://doi.org/10.1016/0364-5916(88)90019-3
- No free full-text PDF was verified. The paper assumes no mutual solid solubility and calculates liquid equilibria with a sub-regular solution model using Miedema interaction parameters.

### Torre et al. (2021)
**Investigation on the Thermodynamic Stability of Nanocrystalline W-Based Alloys: A Combined Theoretical and Experimental Approach.** Materials 14, 7179 (2021).
- Free PDF: https://eprints.soton.ac.uk/454753/1/materials_14_07179_v2.pdf
- DOI: https://doi.org/10.3390/ma14237179
- Reports a strongly positive Ag–W mixing enthalpy of about 40 kJ mol⁻¹ at equiatomic composition, retained phase segregation after prolonged ball milling, and demixing on annealing. This supports immiscibility trends; it is not a complete CALPHAD assessment or a 1:1-mass composite laser study.

### Jiao & Schuh (2018)
**Nanocrystalline Ag-W alloys lose stability upon solute desegregation from grain boundaries.** Acta Materialia 161, 194–206 (2018).
- Free accepted-manuscript PDF already archived: `papers/Jiao_Nanocrystalline_Ag-W_Alloys.pdf`
- DOI: https://doi.org/10.1016/j.actamat.2018.09.014
- Studies only 0.3–1.9 at.% W in nanocrystalline Ag. Useful for grain-boundary segregation and desegregation, not for the bulk 1:1-mass Ag/W phase diagram.

### Tolias (2017) — tungsten thermophysical properties
- Free arXiv PDF: https://arxiv.org/pdf/1703.06302
- Reports (T_m(W)=3695) K and reviews experimental latent heat of fusion clustered around 50 kJ mol⁻¹; 52.3 kJ mol⁻¹ is the commonly cited handbook value. It argues the older 35.3 kJ mol⁻¹ estimate is not supported by measurements.

### NASA high-temperature levitation report (1985)
**Measurement of the Properties of Tungsten at High Temperatures**, NASA TM 87568.
- Free PDF: https://ntrs.nasa.gov/api/citations/19860000586/downloads/19860000586.pdf?attachment=true
- Reports W heat of fusion (52\pm5) kJ mol⁻¹; liquid Ag (C_p=32.6\pm2.1) J mol⁻¹ K⁻¹ over 1281–1549 K and Ag heat of fusion (10.916\pm0.415) kJ mol⁻¹.

### Waseda & Ohtani (1974) — liquid Ag structure
**Static Structure of Liquid Noble and Transition Metals by X-Ray Diffraction.** Physica Status Solidi B (1974).
- DOI: https://doi.org/10.1002/pssb.2220620224
- Wiley PDF endpoint returned access denied; DOI retained.
- Reports X-ray diffraction structure factors for liquid Ag among other liquid metals.

## Access status summary

- Free PDF links verified: Lin 2008, Daraszewicz 2014, Medvedev 2020, Akhmetov 2024, Giri 2023, Byskov thesis, Torre 2021, Tolias 2017, NASA TM 87568. Jiao–Schuh 2018 is already present in the repository.
- DOI-only / no free PDF verified: Vijayakumar et al. 1988 Ag–W phase diagram; Waseda & Ohtani 1974 liquid Ag structure; 2024 Au/Pt–Au/Al electronic interface conductance. Zhang et al. 2004 Ag–W EAM remains DOI-only as recorded in the earlier literature register.
- PDF binary mirroring status: external PDF binaries could not be fetched into this execution environment, so these are verified direct PDF links, not files claimed to have been downloaded into GitHub. Existing user-uploaded PDFs remain archived.
