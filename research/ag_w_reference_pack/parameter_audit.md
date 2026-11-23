# Dai et al. (2006) Ag–W EFS parameter audit

Primary source: X. D. Dai, Y. Kong, J. H. Li, B. X. Liu, “Extended Finnis–Sinclair potential for bcc and fcc metals and alloys,” *Journal of Physics: Condensed Matter* 18 (2006), 4527–4542, DOI [10.1088/0953-8984/18/19/008](https://doi.org/10.1088/0953-8984/18/19/008). The uploaded 17-page PDF is saved as `papers/Extended_Finnis-Sinclair_potential_for_bcc_and_fcc.pdf`.

## Functional form and units

Equations (1)–(5), printed pages 4528–4529, give

```text
U = 1/2 sum_ij V(r_ij) - sum_i f(rho_i)
f(rho_i) = sqrt(rho_i)
rho_i = sum_(j != i) A^2 phi(r_ij)
V(r) = (r-c)^2 * (c0+c1*r+c2*r^2+c3*r^3+c4*r^4), r <= c
phi(r) = (r-d)^2 + B^2*(r-d)^4, r <= d
```

Both radial functions are zero outside their stated cutoff. The parameter tables print A in eV Å⁻¹, d and c in Å, c0 through c4 in eV Å⁻² through eV Å⁻⁶, and B in Å⁻². The JSON transcription previously labeled A as eV¹ᐟ² Å⁻¹. Comparing the uploaded paper corrected that unit to eV Å⁻¹. The numerical parameter values and tabulated potential are unchanged by this correction.

## Transcribed values checked against tables

| Element / pair | A | d (Å) | c (Å) | c0 | c1 | c2 | c3 | c4 | B |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Ag, Table 4, p. 4535 | 0.325514 | 4.41 | 4.76 | 10.6812 | −12.04517 | 5.203072 | −1.013304 | 0.0742308 | −1.293394 |
| W, Table 1, p. 4530 | 1.885948 | 4.41 | 3.25 | 48.52796 | −33.79621 | 5.854334 | −0.0098221 | 0.033338 | 0 |
| Ag–W, Table 6, p. 4540 | 0.6962 | 4.60 | 4.30 | 25.6180 | −26.3407 | 9.1414 | −1.0693 | 0 | 0 |

The LAMMPS `eam/fs` file uses F(rho) = −sqrt(rho), density arrays A² phi(r), and the tabulated pair function r V(r). This follows LAMMPS `eam/fs` setfl conventions and retains the paper’s numerical equations.

## Reference values and implementation checks

Tables 1 and 4 report Ag a = 4.090 Å, Ec = 2.950 eV/atom and W a = 3.160 Å, Ec = 8.916 eV/atom. Table 7 reports AgW B2 and Ag3W L12 benchmarks. Their results and comparison to a generated-table calculation are in `validation/dai_2006_efs_validation_report.json`.

Table 3 reports W melting at 4500 K versus 3695 K experimentally. The paper describes this as the largest melting-point error in its six-metal test, 21.78%. This remains a limitation for high-temperature applications.

The uploaded 17-page Dai article contains no separate supplementary-information citation or file. No Dai supplementary file was included with the upload.
