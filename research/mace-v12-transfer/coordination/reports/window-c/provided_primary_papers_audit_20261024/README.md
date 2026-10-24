# Supplied full-text paper audit — 2026-10-24

## Result

Both supplied PDFs were identified and checked against the candidate structures and source records. Bai et al.'s Table 1 confirms the parameters used to reconstruct the Ti₃SiC₂ CIF. Suh et al.'s paper confirms the room-temperature fcc Ag lattice parameter and provides solid-state thermal expansion data. Its experiment is the same source cited by COD 9013045, so it is not a second independent experimental record.

The PDFs remain in the user-upload attachment area and are not copied into this repository. Their SHA256 values below identify the exact supplied bytes. The Suh file is a third-party scanned PDF; its printed journal pages and tables were legible and visually checked. A duplicate upload of the Suh scan had the same size and SHA256, so it is one file/source.

## Bai et al. (2010): Ti₃SiC₂ DFT structure

Y. Bai, X. He, Y. Sun, C. Zhu, M. Li, and L. Shi, “Chemical bonding and elastic properties of Ti₃AC₂ phases (A = Si, Ge, and Sn): A first-principle study,” *Solid State Sciences* 12 (2010), 1220–1225. DOI: [10.1016/j.solidstatesciences.2010.03.007](https://doi.org/10.1016/j.solidstatesciences.2010.03.007).

The supplied file is the 30-page accepted manuscript. Its Table 1 (printed page 18; PDF page 20) reports Ti₃SiC₂ as space group P6₃/mmc, `a = 3.0592 Å`, `c = 17.636 Å`, `zTi2 = 0.1351`, and `zC = 0.5726`. Visual inspection confirms the table columns and values. The candidate CIF uses these exact values and the corresponding Wyckoff sites; the CIF is a reconstructed derivative, not a source-original file. The methods identify CASTEP plane-wave pseudopotential DFT with GGA-PW91 and full optimization of lattice parameters and internal coordinates.

This verifies the numerical link between the supplied CIF and an independent computational paper. It does not make the geometry a distinct polymorph or an independent experimental record. The six derived supercell files remain one geometry-only source family; no training, development, or test role is assigned.

## Suh, Ohta & Waseda (1988): Ag thermal expansion

I.-K. Suh, H. Ohta, and Y. Waseda, “High-temperature thermal expansion of six metallic elements measured by dilatation method and X-ray diffraction,” *Journal of Materials Science* 23 (1988), 757–760. DOI: [10.1007/BF01174717](https://doi.org/10.1007/BF01174717).

The supplied four-page scan includes the complete article. In Table II (printed page 758), silver is fcc and has an X-ray lattice parameter of `0.4086 nm` at 293 K (`4.086 Å`), matching COD 9013045. The same table reports `0.4172 nm` at 1228 K with `ΔL/L = 2.13%`; a second run reports `0.4083 nm` at 293 K and `0.4175 nm` at 1226 K. The authors state that Ag remains face-centred cubic over the quoted temperature range.

Table III (printed page 759) gives the X-ray fit over 293–1228 K:

`ΔL/L = A1(T − 293) + A2(T − 293)²`, with `A1 = 1.81 × 10⁻⁵ K⁻¹` and `A2 = 4.74 × 10⁻⁹ K⁻²`.

COD 9013045 cites this same 1988 experiment as its source. The article directly verifies the experimental provenance behind the COD record; the COD entry and paper must remain one experimental lineage in any source count or validation split. These are solid-Ag thermal expansion data and do not provide a liquid-Ag density anchor, so the existing liquid density reference and liquid candidate volumes are unchanged.

## Integrity and scope

| Supplied file | Pages | SHA256 |
|---|---:|---|
| `AcceptedManuscript.pdf` | 30 | `cf393d8783f8d87da4b35ed7093d494b6deb38b20edc4bb8eea81d94efd2f755` |
| `High-temperature thermal expansion of six metallic elements measured by dilatation method and X-ray diffraction(科研通-ablesci.com).pdf` | 4 | `f4dc3c72a5ec03520a53815a1d44f0e5c552872d363946ca34760cd06ecc1911` |

Checks were limited to bibliographic identity, relevant tables, and their relation to the existing candidate geometries. No DFT, molecular dynamics, model training, density-anchor change, or dataset-role assignment was performed.
