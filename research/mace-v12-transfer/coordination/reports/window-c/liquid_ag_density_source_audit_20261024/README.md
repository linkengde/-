# Liquid Ag density source search — 2026-10-24

## Candidate references found

**Abdullaev et al. (2020):** R. N. Abdullaev, R. A. Khairulin, and S. V. Stankus, “Density and thermal expansion of silver in the solid and liquid states,” *Journal of Physics: Conference Series* 1677 (2020), 012161. DOI: [10.1088/1742-6596/1677/1/012161](https://doi.org/10.1088/1742-6596/1677/1/012161). Crossref and OpenAlex metadata agree on the title, authors, journal, year, volume, and article number. OpenAlex marks the published version as open access under CC BY. Its abstract says pure silver density and thermal expansion were investigated from 293 to 1631 K, including measurements through the solid–liquid transition, and that volumetric-property errors were estimated.

**Korobeinikov et al. (2018):** I. I. Korobeinikov, D. Chebykin, X. Yu, S. Seetharaman, and O. Volkova, “Density of tin, silver and copper,” *Archives of Materials Science and Engineering* 92(1) (2018), 28–32. DOI: [10.5604/01.3001.0012.5509](https://doi.org/10.5604/01.3001.0012.5509). Crossref and OpenAlex identify the same article. The abstract describes Archimedean measurements of liquid-metal densities from 700 to 1520°C, a range that includes 1250 K (976.85°C). It reports deviations of 0.5–1% from literature values depending on the metal; that is an agreement statement, not a pointwise measurement uncertainty. This is a promising source for the 1250 K reference if the article's silver table includes a nearby point and reports usable uncertainty.

## Current audit status

**Status: promising sources identified; the liquid-density gate is still open.** The IOP article page/PDF and the ArchivesMSE PDF returned HTTP tunnel error 403 from this environment. No silver table, equation, 1250 K interpolation, numerical density, or pointwise uncertainty was extracted. The existing model-derived density anchor remains unchanged, and no liquid-Ag candidate volumes were rescaled.

The cloud-environment configuration draft now includes `iopscience.iop.org` and `archivesmse.org`, but the draft is not active until it is reviewed, saved, and published in environment settings. Once active, retry the publisher PDFs and record the exact silver data row, 1250 K value or interpolation method, and uncertainty before changing candidate volumes.

## Reproduction notes

- Crossref API record: `https://api.crossref.org/works/10.1088%2F1742-6596%2F1677%2F1%2F012161`
- OpenAlex API record: `https://api.openalex.org/works/https://doi.org/10.1088/1742-6596/1677/1/012161`
- Publisher PDF URL: `https://iopscience.iop.org/article/10.1088/1742-6596/1677/1/012161/pdf`
- Crossref API record: `https://api.crossref.org/works/10.5604%2F01.3001.0012.5509`
- OpenAlex API record: `https://api.openalex.org/works/https://doi.org/10.5604/01.3001.0012.5509`
- Publisher PDF URL: `https://archivesmse.org/gicid/pdf/01.3001.0012.5509`
- Search results alone are not accepted as a numerical reference; no value has been inferred from the abstract.
