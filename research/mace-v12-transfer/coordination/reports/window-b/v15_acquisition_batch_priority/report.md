# Initial future acquisition batch — Window B

Recommend six existing label-free proposals: two per interface, one wider distance endpoint and one passing registry shift. A decides DFT after current screening. Geometry PASS is triage, not chemical stability or independent morphology validation. No candidates were regenerated or changed; no DFT, model inference/training, MD/TTM or V15 blind-label access occurred.

| Priority | Candidate | Axis/offset Å | Pair Å | Backup |
|---|---|---|---|---|
| 1 | `AgC_wide_distance_p0p35` | distance / +0.35 | 2.6100 | AgC_wide_distance_m0p20 |
| 2 | `AgC_wide_registry_p0p45` | registry / +0.45 | 2.3043 | No second passing registry |
| 3 | `AgSi_wide_distance_p0p50` | distance / +0.50 | 2.9401 | AgSi_wide_distance_m0p20 |
| 4 | `AgSi_wide_registry_m0p30` | registry / -0.30 | 2.4585 | AgSi_wide_registry_p0p20 |
| 5 | `AgTi_wide_distance_p0p35` | distance / +0.35 | 2.9380 | AgTi_wide_distance_m0p20 |
| 6 | `AgTi_wide_registry_p0p45` | registry / +0.45 | 2.6268 | No second passing registry |

All three distance selections expand the historical same-motif PW upper-distance range; registry selections add the lateral axis rather than selecting another nearby normal point. AgSi needs wider surrounding coverage (V14 local/outside RMSE 0.1266/0.1079 eV/Å), but all existing shell variants are near their training parents; defer them from this first batch. AgC/AgTi errors concentrate near the contact shell. V14 force and separation errors guide geometry axes only; no unseen V15 scores or blind outcomes influenced selection.

Each JSON selection preserves exact file/geometry hashes, parent frame/file hash/marked IDs, translation, novelty and all element-pair minima. Backup distance points are alternate compression samples; they are not paired with the selected endpoint by default. AgSi has a second passing registry (+0.20) as backup; AgC/AgTi have only one passing registry each, so no flagged registry is promoted.

## All candidate near-pairs

| Pair | Decision |
|---|---|
| `AgC_wide_distance_m0p35` / `AgC_wide_distance_m0p20` | Neither selected; defer flagged contact/near-parent shell variants. |
| `AgC_wide_distance_p0p20` / `AgC_wide_distance_p0p35` | Keep selected endpoint only; defer its adjacent near partner. |
| `AgC_wide_surroundings_m0p12` / `AgC_wide_surroundings_p0p12` | Neither selected; defer flagged contact/near-parent shell variants. |
| `AgSi_wide_distance_m0p35` / `AgSi_wide_distance_m0p20` | Neither selected; defer flagged contact/near-parent shell variants. |
| `AgSi_wide_distance_p0p20` / `AgSi_wide_distance_p0p35` | Neither selected; defer flagged contact/near-parent shell variants. |
| `AgSi_wide_registry_m0p45` / `AgSi_wide_registry_m0p30` | Keep selected endpoint only; defer its adjacent near partner. |
| `AgSi_wide_distance_p0p35` / `AgSi_wide_distance_p0p50` | Keep selected endpoint only; defer its adjacent near partner. |
| `AgSi_wide_surroundings_m0p12` / `AgSi_wide_surroundings_p0p12` | Neither selected; defer flagged contact/near-parent shell variants. |
| `AgTi_wide_distance_m0p35` / `AgTi_wide_distance_m0p20` | Neither selected; defer flagged contact/near-parent shell variants. |
| `AgTi_wide_distance_p0p20` / `AgTi_wide_distance_p0p35` | Keep selected endpoint only; defer its adjacent near partner. |
| `AgTi_wide_surroundings_m0p12` / `AgTi_wide_surroundings_p0p12` | Neither selected; defer flagged contact/near-parent shell variants. |

Reviewed all 11 listed near pairs: zero pairs have both candidates in the initial six. Defer adjacent +0.20/+0.35 C/Ti distance samples and +0.35 Si in favor of the wider selected endpoint. Si distance +0.20/+0.35 is a chain; choose +0.50 first, not the complete chain. Short flagged compression/registry pairs remain flagged. Opposite ±0.12 shell variants remain near-parent acquisitions, not independent tests.

## Evidence limits and A handoff

All 33 candidate byte hashes match the design manifest; six selections are label-free with preserved original screening PASS. The 14 passing proposals have existing-reference exact/near screen counts zero in the original 55-file/166-frame inventory, including residual acquisitions and frozen V15 INPUTS. This report reuses that inventory and does not claim a fresh all-repository screen. A must rescreen selected/backup geometries against any subsequently added train, test, planned acquisition or frozen blind INPUT before assigning DFT; preserve blind labels unused.

V15 provenance remains limited: 14 partial records and 12 explicitly declared inherited rows lack complete original-run revalidation. Future labels should use matched, recorded calculator parameters and energy convention; original-history recovery/relabeling is a separate job. Do not treat passing pair distances, rank or acquisition selection as model reliability.

Inputs/hashes are in batch_recommendation.json and SHA256SUMS.txt. Structures are referenced at existing paths and not duplicated.
