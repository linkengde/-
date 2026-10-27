# Jeitschko–Nowotny 1967 Ti₃SiC₂ source audit

## Result

The uploaded full paper confirms an independent single-crystal experimental source for Ti₃SiC₂: W. Jeitschko and H. Nowotny, *Monatshefte für Chemie* 98 (1967), 329–337, DOI `10.1007/BF00899949`. It is distinct from the 1998 Kisi et al. neutron-diffraction source associated with COD 9009647 / ICSD 86213. The paper provides enough crystallographic information to construct and screen a derivative CIF; an official ICSD 25762 export is not required to establish this paper’s source lineage.

The supplied package's CIF is a faithful reconstruction from the rounded abstract values (`zTi=0.135`, `zC=0.5675`). The newly supplied full paper also reports later Fourier estimates (`zTi=0.1357(15)`, `zC=0.570(4)`). Both variants are retained here as separate candidates; neither is represented as an original ICSD export. The abstract-value CIF expands to 12 atoms, has cell volume `144.03011684 Å³`, and a minimum periodic distance of `2.135409 Å`. The full-text Fourier variant has the same cell and composition, and minimum periodic distance `2.117811 Å`.

The text also discusses possible Ti-site deficiency / Ti–Si substitution and likely carbon vacancies. Consequently, unit occupancies in these idealized CIFs are a modeling assumption, not a resolved experimental occupancy result.

## Evidence and file identity

- Uploaded original-paper PDF SHA256: `6765f99b3cc83daf3e48e4674dbeba109e5fffef8fa66cd53e266733ac4a4a63`.
- Uploaded ZIP SHA256: `09c311c3a1af2845031e9bf0978030b27684ec67f8493d3e86de405f0ac2c10b`.
- All nine entries in the ZIP's `SHA256SUMS.txt` passed verification.
- The uploaded abstract-value CIF is preserved byte-for-byte. Its SHA256 is `5fd657f81d55bb0d15e4c2e12911846f537d97f28244f53952eee23e6d1edb5f`.
- The paper's page 333 lists the four Wyckoff site sets; pages 334–335 report the refined Fourier estimates and discuss occupancy uncertainty. The original paper does not itself identify an ICSD accession. The reported mapping to ICSD 25762 remains secondary-source metadata until an official record is obtained.

## Dataset role

Register the source as `source_candidate / independent_experiment / role_pending`. This establishes an independent publication/source family, but does not itself establish an independent development or test entry. Before assigning a data role, screen the actual derived interface geometries against the complete training candidate and preserve parent/source lineage. A single crystallographic parent does not supply three independent parent families across Ag–C, Ag–Si, and Ag–Ti contacts.

The CIFs provide geometry only. They contain no DFT energies or forces and do not resolve the failed V19 development gates. Do not use them as training/dev labels. DFT execution and model fitting remain blocked by the existing resource/method review and registered-role requirements.

## Reproduction

Run `python verify_structure.py` with ASE installed to parse both CIFs and check formula, atom count, space group, cell, volume, and minimum periodic distance. This is a geometry-only check; it does not validate occupancies, source-database identity, or model suitability.
