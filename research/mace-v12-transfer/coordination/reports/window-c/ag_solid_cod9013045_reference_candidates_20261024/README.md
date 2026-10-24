# Solid Ag source and geometry candidates — 2026-10-24

The previous B1 solid-Ag proposal used `a = 4.09 Å` from the legacy Ag EAM model. This package adds an experimentally sourced geometry family using the 293 K fcc Ag record in COD 9013045: `a = b = c = 4.086 Å`, space group Fm-3m (No. 225), DOI [10.1007/BF01174717](https://doi.org/10.1007/BF01174717). The COD record cites I.-K. Suh, H. Ohta, and Y. Waseda, “High-temperature thermal expansion of six metallic elements measured by dilatation method and X-ray diffraction,” *Journal of Materials Science* 23 (1988), 757–760. The user supplied the complete four-page article scan. Table II directly confirms the same 293 K Ag value (`0.4086 nm = 4.086 Å`), the fcc structure, and solid-state Ag lattice expansion; Table III gives the X-ray fit through 1228 K. This verifies the paper behind the COD record, not a second independent experiment. The COD CIF was retrieved from its official record (4,769 bytes; SHA256 `6dbb7ae83a23963204724bb2da599fbe6304a6c875d5cc3f28429639a582fe31`). Its header credits the American Mineralogist Crystal Structure Database and requires attribution; the source CIF itself is not redistributed here. PDF page/table details and the uploaded PDF hash are in [the supplied-paper audit](../provided_primary_papers_audit_20261024/README.md).

`candidates/` contains six reproducible Ag geometries at the COD lattice parameter: a 32-atom conventional 2×2×2 cell and a 64-atom primitive 4×4×4 cell, each with its parent plus two seeded 0.025 Å Gaussian displacement variants. All six are one source family and can be considered only as training candidates after B role review. They are not independent polymorphs, validation structures, or DFT labels. The former 4.09 Å EAM-derived candidates remain unchanged.

The generator performs extxyz read-back, formula, source hash, finite geometry, PBC, minimum-distance, unique-file-hash, and DFT-disabled checks. No calculation was launched and no dataset or task card was modified.

Reproduce from the repository root:

```bash
/workspace/.venvs/mace/bin/python research/mace-v12-transfer/coordination/reports/window-c/ag_solid_cod9013045_reference_candidates_20261024/build_candidates.py
```
