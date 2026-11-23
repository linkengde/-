# XTANT-3 Ag/W TTM tables — provenance and use limits

Copied from the public author repository `N-Medvedev/XTANT-3_coupling_data`, ref `TTM_parameters_v3`, on 2026-11-23. Original text columns and units are retained.

## Files and grids

| Material | Quantity | Source file | Grid | Units |
|---|---|---|---:|---|
| W | electron heat capacity (C_e(T_e)) | `Ce_metals/Ce_W.txt` | 500–19,000 K | J m⁻³ K⁻¹ |
| W | electron-ion coupling (G(T_e)) | `G_e-ph_metals/G_e-ph_W.txt` | 500–19,000 K | W m⁻³ K⁻¹ |
| W | electronic thermal conductivity | `K_metals/K_W.dat` | 300–50,000 K | W m⁻¹ K⁻¹ |
| Ag | electron heat capacity (C_e(T_e)) | `Ce_metals/Ce_Ag.txt` | 600–23,000 K | J m⁻³ K⁻¹ |
| Ag | electron-ion coupling (G(T_e)) | `G_e-ph_metals/G_e-ph_Ag.txt` | 600–23,000 K | W m⁻³ K⁻¹ |
| Ag | electronic thermal conductivity | `K_metals/K_Ag.dat` | 300–50,000 K | W m⁻¹ K⁻¹ |

The (K) files contain columns `Te, kappa_tot, kappa_e_ph, kappa_e_e, mu, Ce`. The total conductivity is assembled from the electron-phonon and electron-electron channels using Matthiessen's rule. The (C_e) column inside each (K) file is **not the same calculation** as the separate `Ce_*.txt` file: the separate (C_e) data are from dynamical simulations at the gamma point and a given atomic temperature, while the (K)-file (C_e) is for an ideal crystal averaged over a 7x7x7 k-point mesh. Do not silently mix these two (C_e) definitions.

## Important scope limits

1. The dataset is intended for **high electronic temperatures**. Its Zenodo description explicitly warns that room-temperature values may not always be reliable.
2. (G) is tabulated as a function of (T_e) only. It does not provide a complete (G(T_e,T_l)) surface.
3. The (K) files are tabulated versus (T_e), not a two-dimensional (k_e(T_e,T_l)) surface. They are useful high-(T_e) transport references, but do not by themselves close the requested (T_e,T_l)-dependent conductivity.
4. The W values are not a single uncontested truth. For example, the XTANT-3 W (G) table starts at (1.2711\times10^{16}) W m⁻³ K⁻¹ at 500 K, whereas Daraszewicz et al. (2014) report an ab-initio (G_0=1.65\times10^{17}) and a pump-probe fit (1.4(3)\times10^{17}) W m⁻³ K⁻¹. Lin et al. (2008) also report substantially different low-(T_e) estimates and a rising (G(T_e)). Preserve these as separate model/data branches and run sensitivity comparisons; do not average them into one curve without resolving the method differences.
5. For W (k_e(T_e,T_l)), Lin et al. (2008) describe a Drude form (k_e=v_F^2 C_e(T_e)\tau_e(T_e,T_l)/3), with (1/\tau_e=A T_e^2+B T_l). Daraszewicz et al. (2014) use a practical Drude parameterization (k_e(T_e,T_l)=k_{e0}[C_e(T_e)/C_e(T_0)](T_0/T_l)^n), (n=1.5), (T_0=300) K. Their 300 K total W conductivity is 174 W m⁻¹ K⁻¹ and (k_e/k=0.75), implying (k_{e0}=130.5) W m⁻¹ K⁻¹ for that model. This is a model closure, not a directly measured full 2D surface.
6. These are elemental Ag and W tables. They do not define mixed-cell Ag/W electronic transport, interface electron conductance, or the partition of electronic heat capacity in a mixed finite-volume cell.

## Source provenance

- Author repository: https://github.com/N-Medvedev/XTANT-3_coupling_data
- Dataset ref: `TTM_parameters_v3`
- Zenodo record: https://zenodo.org/records/10183636
- Zenodo DOI shown on the record: `10.5281/zenodo.10183636`
- The repository's `CITATION.cff` lists `10.5281/zenodo.10183635`; retain this citation discrepancy in provenance records rather than silently rewriting it.

Source Git blob SHAs:
- `Ce_metals/Ce_W.txt`: 75b069d04fdb0c77f162b522e6c1ebae3e0b5a8f
- `G_e-ph_metals/G_e-ph_W.txt`: 50a7eaaceb89b4243dfcffb88dadbca3c066f5bd
- `K_metals/K_W.dat`: 15776fd49250b9d87488dbcbcf272f043fdf7867
- `Ce_metals/Ce_Ag.txt`: a7bf7724531bae6e75298d1da480f29cd3f37272
- `G_e-ph_metals/G_e-ph_Ag.txt`: 030d7cf6d2797055e2cd0a1e980935478102e175
- `K_metals/K_Ag.dat`: 3ad8513910533f0f1925badfe663aa3480211152

No TTM or MD trajectory was run as part of this data import.
