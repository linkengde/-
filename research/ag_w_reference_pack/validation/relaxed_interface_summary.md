# Relaxed Ag/W interface screening with the Dai EFS potential

## Scope and method

This is a preliminary zero-temperature LAMMPS screening of Ag(100)/W(100) and Ag(111)/W(110) using the generated `Dai_2006_EFS_AgW.eam.fs` potential (SHA-256 `215c53d03029736f3be4de5f8833a0d5a0e99fd18cd72b9ba118a0b583ff04f0`). It is not a converged reproduction of the DFT calculations in He et al. (2018), Table 4.

Each periodic interface cell contains two interfaces. The Ag(100)/W(100) cell contains 64 Ag and 112 W atoms; the Ag(111)/W(110) cell contains 1080 Ag and 1008 W atoms. The in-plane dimensions and cell height are held fixed at each scan point. Central bulk-like slab cores are fixed, and the other atoms are minimized with conjugate gradients (energy tolerance 0, force tolerance 1e-5 eV/Å). Initial gaps of 2.0, 2.4, 2.8, 3.2, 3.6, and 6.0 Å were sampled for each lateral registry. The observed maximum residual force at the selected minimum is below 1.81e-6 eV/Å. Bulk reference energies use the same in-plane strain and a relaxed normal repeat length; surface references use the same strained slabs and constraints.

The reported quantities are

\[
W_{sep}=\frac{E_{Ag,surf}+E_{W,surf}-E_{interface}}{2A},\qquad
\gamma_{interface}=\frac{E_{interface}-N_{Ag}e_{Ag,bulk}-N_We_{W,bulk}}{2A}.
\]

The separated 6 Å scan agrees with the sum of the two independently minimized surface slabs within 7e-10 eV. The identity \(W_{sep}=\gamma_{Ag,surf}+\gamma_{W,surf}-\gamma_{interface}\) closes to within 7e-12 J/m².

## Results

| Interface | Registry | Lowest-energy initial gap (Å) | (W_{sep}) (J/m²) | γ interface (J/m²) | He et al. DFT (W_{sep}) / γ (J/m²) |
|---|---:|---:|---:|---:|---:|
| Ag(100)/W(100) | top | 2.8 | 2.894 | 1.127 | 2.80 / 0.72 |
| Ag(100)/W(100) | bridge x | 2.4 | 1.654 | 2.368 | 2.80 / 0.72 |
| Ag(100)/W(100) | bridge diagonal | 2.0 | 2.609 | 1.412 | 2.80 / 0.72 |
| Ag(111)/W(110) | origin | 2.8 | 3.192 | 0.156 | 2.45 / 0.86 |
| Ag(111)/W(110) | shift x | 2.8 | 3.214 | 0.134 | 2.45 / 0.86 |
| Ag(111)/W(110) | shift xy | 2.8 | 3.206 | 0.142 | 2.45 / 0.86 |
| Ag(111)/W(110) | shift mixed | 2.8 | 3.222 | 0.127 | 2.45 / 0.86 |

The lowest-energy sampled Ag(100)/W(100) registry gives a work of separation 0.094 J/m² above the DFT value and an interface energy 0.407 J/m² above it. The best sampled Ag(111)/W(110) registry gives a work of separation 0.772 J/m² above DFT and an interface energy 0.733 J/m² below DFT. The second orientation therefore remains a material disagreement for this EFS model. The registry spread is large for Ag(100)/W(100) and small for the sampled Ag(111)/W(110) cell.

The matched cell strains are substantial for Ag(100)/W(100): Ag +4.53% and W -4.33% in both lateral directions. For Ag(111)/W(110), the strains are Ag +2.02% and -0.83%, and W -1.98% and +0.84%. The potential surface energies are 0.910/3.112 J/m² (Ag/W, 100/100) and 0.730/2.618 J/m² (Ag/W, 111/110). Strain, termination, core constraints, and the exact match cell can all affect comparison with DFT.

## Interpretation

The Ag(100)/W(100) top registry is a useful positive check against He et al., but one close value does not validate the full potential. The Ag(111)/W(110) result and the strong Ag(100)/W(100) registry dependence show that the model still needs more interface and thermodynamic checks before production MD. This calculation does not include finite-temperature dynamics and does not authorize high-temperature transport or mechanical-property claims.

The full gap-by-registry data, energies, strains, constraints, residual forces, reference checks, LAMMPS version, and potential hash are in `relaxed_interface_report.json`.
