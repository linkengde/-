# User-supplied XTANT-3 Ag/W tables

This folder keeps the original Ag/W-only `Ce`, `G_e-ph`, and `K` text tables supplied for the Ag–W TTM review. The data files are copied without editing. The uploaded upstream README is preserved as `SOURCE_README.md`; `LICENSE` and `.gitattributes` are preserved from the same upload batch. SHA-256 hashes are listed in `../SHA256SUMS.txt`.

| File | Role | Supplied temperature range |
|---|---|---:|
| `Ce_Ag.txt`, `Ce_W.txt` | Dynamical-simulation electron heat capacity | Ag 600–23,000 K; W 500–19,000 K |
| `G_e-ph_Ag.txt`, `G_e-ph_W.txt` | Electron-phonon/electron-ion coupling as a function of `Te` | Ag 600–23,000 K; W 500–19,000 K |
| `K_Ag.dat`, `K_W.dat` | Total electronic thermal conductivity, channel contributions, chemical potential, and a separate `Ce` value | Both 300–50,000 K |

Do not substitute the `Ce` columns inside `K_Ag.dat` and `K_W.dat` for the standalone `Ce_*.txt` files without documenting that choice. According to `SOURCE_README.md`, these use different calculations: ideal crystal with a 7×7×7 k-point mesh for the K-table `Ce`, versus dynamical simulations at the gamma point for standalone `Ce`.

The source README describes the data as GPL-3.0 and also states “non-commercial peaceful purposes only.” Both source notices are retained here. This folder does not provide a legal interpretation; preserve the notices and confirm applicable terms before redistribution or commercial use.

The full comparison, temperature-grid inventory, and TTM limits are in [`../XTANT3_AG_W_DATA_AUDIT_20261124.md`](../XTANT3_AG_W_DATA_AUDIT_20261124.md). No files for unrelated metals or semiconductors are included in this Ag–W-specific subfolder.
