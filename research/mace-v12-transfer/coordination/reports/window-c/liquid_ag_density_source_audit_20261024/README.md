# Liquid Ag density source search — 2026-10-24

## Candidate references found

**Abdullaev et al. (2020):** R. N. Abdullaev, R. A. Khairulin, and S. V. Stankus, “Density and thermal expansion of silver in the solid and liquid states,” *Journal of Physics: Conference Series* 1677 (2020), 012161. DOI: [10.1088/1742-6596/1677/1/012161](https://doi.org/10.1088/1742-6596/1677/1/012161). Crossref and OpenAlex metadata agree on the title, authors, journal, year, volume, and article number. OpenAlex marks the published version as open access under CC BY. Its abstract says pure silver density and thermal expansion were investigated from 293 to 1631 K, including measurements through the solid–liquid transition, and that volumetric-property errors were estimated.

**Korobeinikov et al. (2018):** I. I. Korobeinikov, D. Chebykin, X. Yu, S. Seetharaman, and O. Volkova, “Density of tin, silver and copper,” *Archives of Materials Science and Engineering* 92(1) (2018), 28–32. DOI: [10.5604/01.3001.0012.5509](https://doi.org/10.5604/01.3001.0012.5509). Crossref and OpenAlex identify the same article. The abstract describes Archimedean measurements of liquid-metal densities from 700 to 1520°C, a range that includes 1250 K (976.85°C). It reports deviations of 0.5–1% from literature values depending on the metal; that is an agreement statement, not a pointwise measurement uncertainty. OpenAlex marks the article closed-access, so the silver table and uncertainty still need full-text inspection.

**Assael et al. (2012):** “Reference data for the density and viscosity of liquid antimony, bismuth, lead, nickel and silver,” *High Temperatures–High Pressures*. OpenAlex record [W2608312387](https://openalex.org/W2608312387) lists M. J. Assael, Jiangtao Wu, Erhard Kaschnitz, Konstantinos D. Antoniadis, W. A. Wakeham, Agni E. Kalyva, R. Michael Banish, and I. Egry, and points to the Aristotle University of Thessaloniki institutional repository record at `https://ikee.lib.auth.gr/record/217743`. OpenAlex supplies no DOI or PDF and marks the record closed-access. An earlier exact-title Crossref search did not resolve a DOI; a follow-up request on 2026-10-24 returned HTTP 429, so DOI resolution remains open.

### Assael source follow-up (2026-10-24)

- The user reports that the [publisher's official issue page](https://www.oldcitypublishing.com/journals/hthp-home/hthp-issue-contents/hthp-volume-41-number-3-2012/hthp-41-3-p-161-184/) confirms the complete title, authors, *High Temperatures–High Pressures* 41(3), 161–184, and provides a “full text (IP)” entry. The page abstract reports a 0.9% standard deviation for the recommended reference correlation at 95% confidence. This is not treated as a 1250 K single-point measurement uncertainty. The publisher page request from this session was blocked by the network tunnel (HTTP 403), so these page details are recorded as user-verified, not independently fetched here.
- The user reports that the [OpenArchives aggregator record](https://www.openarchives.gr/aggregator-openarchives/edm/IKEE_AUT/000097-217743?language=en) marks the IKEE item `openAccess`, while its “view original record” link currently returns a gateway error. This is an access-status conflict between an aggregator flag and the repository/publisher evidence; it does not establish that a full-text file is retrievable. This session's OpenArchives request was blocked by the network tunnel (HTTP 403), so the aggregator page itself has not been independently inspected.
- The IKEE [original repository record](https://ikee.lib.auth.gr/record/217743?ln=en) was fetched successfully (HTTP 200). It confirms the title, authors, 2012 publication date, and journal citation 41(3), 161–184. The record is labeled “Batch imported”; its “Files” tab is disabled and the file panel is empty. No attached manuscript or downloadable file link is present in the returned record.
- The publisher's “full text (IP)” route and the OpenArchives original-record link therefore remain unverified from this environment. No paper full text, silver correlation equation, correlation temperature interval, 1250 K value/interpolation, or pointwise uncertainty has been obtained. The existing model-derived density anchor remains unchanged, and no liquid-Ag candidate volumes were rescaled.

## Current audit status

**Status: promising citations identified; the liquid-density gate is still open.** A fresh request to the IOP article page and PDF returned HTTP tunnel error 403. The ArchivesMSE PDF endpoint returned HTTP 200 with `text/html` rather than a PDF; its HTML title was “Publishers Panel.” OpenAlex marks that article closed-access. For Assael et al. 2012, the IKEE landing page is accessible but contains metadata only, with the file tab disabled; OpenAlex marks the work closed and provides no PDF. The user-reported OpenArchives `openAccess` flag remains unverified and conflicts with the inaccessible original-record link. No silver table, equation, 1250 K interpolation, numerical density, or pointwise uncertainty was extracted. The existing model-derived density anchor remains unchanged, and no liquid-Ag candidate volumes were rescaled.

The cloud-environment configuration draft includes `iopscience.iop.org`, `archivesmse.org`, `ikee.lib.auth.gr`, `oldcitypublishing.com`, `www.oldcitypublishing.com`, `openarchives.gr`, and `www.openarchives.gr`; these domains are not active until the draft is reviewed, saved, and published in environment settings. Once active, retry the publisher's IP full-text entry and the OpenArchives original-record path, and record the exact silver correlation, validity interval, 1250 K value or interpolation method, and applicable uncertainty before changing candidate volumes.

The IKEE repository record has now been inspected and contains no file. After configuration publication, retry the publisher and aggregator links; do not infer full-text availability from OpenAlex or an aggregator's `openAccess` label alone.

## Reproduction notes

- Crossref API record: `https://api.crossref.org/works/10.1088%2F1742-6596%2F1677%2F1%2F012161`
- OpenAlex API record: `https://api.openalex.org/works/https://doi.org/10.1088/1742-6596/1677/1/012161`
- Publisher PDF URL (latest request: tunnel HTTP 403): `https://iopscience.iop.org/article/10.1088/1742-6596/1677/1/012161/pdf`
- Crossref API record: `https://api.crossref.org/works/10.5604%2F01.3001.0012.5509`
- OpenAlex API record: `https://api.openalex.org/works/https://doi.org/10.5604/01.3001.0012.5509`
- Publisher PDF URL (latest request: HTTP 200 `text/html`, title “Publishers Panel”; OpenAlex says closed access): `https://archivesmse.org/gicid/pdf/01.3001.0012.5509`
- OpenAlex API record for Assael et al. 2012: `https://api.openalex.org/works/W2608312387`
- OpenAlex public record for Assael et al. 2012: `https://openalex.org/W2608312387`
- Publisher issue page with “full text (IP)” entry (user-confirmed page details; current-session request blocked by network tunnel): `https://www.oldcitypublishing.com/journals/hthp-home/hthp-issue-contents/hthp-volume-41-number-3-2012/hthp-41-3-p-161-184/`
- OpenArchives aggregator record (user reports `openAccess`; current-session request blocked by network tunnel): `https://www.openarchives.gr/aggregator-openarchives/edm/IKEE_AUT/000097-217743?language=en`
- IKEE institutional-repository metadata record (fetched HTTP 200; no file attached): `https://ikee.lib.auth.gr/record/217743?ln=en`
- Crossref exact-title query (follow-up request returned HTTP 429): `https://api.crossref.org/works?query.title=Reference%20data%20for%20the%20density%20and%20viscosity%20of%20liquid%20antimony%2C%20bismuth%2C%20lead%2C%20nickel%20and%20silver`
- Search results alone are not accepted as a numerical reference; no value has been inferred from the abstract.
