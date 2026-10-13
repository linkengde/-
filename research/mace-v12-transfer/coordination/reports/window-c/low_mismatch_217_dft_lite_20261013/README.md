# Ag–Ti₃SiC₂ 217-atom low-mismatch DFT-lite structure

This directory publishes the unchanged, geometry-only 217-atom candidate for a future fixed-geometry DFT single point.

- Composition: Ag84Ti57Si38C38
- Atom count: 217
- Ag/Ti₃SiC₂ twist: 12.51983°
- Ti₃SiC₂ in-plane lattice: retained from the audited low-mismatch parent
- Uniform Ag in-plane strain: +1.1029%; Ti₃SiC₂ strain: 0%
- Interface: Ag(111) over Si-terminated Ti₃SiC₂(0001)
- Geometry: unrelaxed; original registry and coordinates are retained, with only layer selection and a uniform z translation/cell shortening to preserve a 20 Å vacuum
- Cell PBC: `[true, true, true]`
- Intended use: fixed-geometry PBE/PAW DFT single-point validation; this package is not a DFT result or convergence claim

Structure SHA-256: `945c5e214c1ceaeba4125ff785f4092e14858c75a985061a7841da3e5ce6a9a2`

The geometry metrics were copied from the generating workspace. `SHA256SUMS.txt` covers the files in this published directory. The previous execution preflight held the DFT because the active 32 GiB cgroup did not have adequate safe memory headroom. No DFT, MD, relaxation, or training was run for this publication. Do not treat the geometry candidate as an energy/force label.
