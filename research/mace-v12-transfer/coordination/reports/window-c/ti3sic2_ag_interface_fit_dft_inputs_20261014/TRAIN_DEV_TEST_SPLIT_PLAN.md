# Proposed TRAIN / DEV / TEST split

Status: proposal only. No DFT labels were generated and no existing data role was changed.

| Parent family | Included structures | Proposed assignment | Restriction |
|---|---|---|---|
| `TSC_KISI_COD9009647_BASELINE_DERIVED_FAMILY` | 4 volume/shear candidates and 2 Si ±0.02 Å candidates; ancestry also includes the baseline and E5 Ti/C response labels | `TRAIN_CANDIDATE` as one indivisible family, pending explicit role freeze | The five E5 numerical-control labels remain numerical-control/diagnostic. They are not copied or silently converted. Because this family was used in earlier diagnostics, none of its descendants can be DEV or TEST. |
| `AG_TSC_SI_TERM_TWIST12P51983_217_PARENT_FAMILY` | 9 rigid Ag–slab gap-scan geometries | `TRAIN_CANDIDATE` for a future unified Ag–Ti–Si–C model only | All nine gaps stay in one family and one split. The family is not independent of its 217-atom parent or the 394-atom parent SHA recorded in its metadata. Not a TSCBOP target. |
| Independent Ti₃SiC₂ parent | No qualifying geometry identified by E7 | `DEV: SOURCE_PENDING` | Must be a genuinely separate source/parent, not another paper's measurement of the same prototype or a perturbation/supercell of Kisi/COD. |
| Independent Ag/Ti₃SiC₂ parent | No complete independently sourced coordinate set with reproducible registry found in the reviewed package | `DEV: SOURCE_PENDING` | Zhang et al. interface energetics can motivate construction but the archived evidence did not contain the full coordinates/forces needed to label this family. |
| Evaluator-held independent parents | Not inspected or enumerated | `TEST: EVALUATOR_CONTROLLED` | C must not access sealed coordinates or labels. Evaluator must certify ancestry disjoint from all training and development families. |

## Split policy

Split only by parent family and all descendants. Never random-split gap points, perturbations, strain signs, k-point variants, surface terminations derived from a common parent, or trajectory frames. A future combined model may train on both candidate TRAIN families, but the three-element TSCBOP fit may consume only Ti–Si–C labels; Ag-containing interface data cannot constrain its parameters. DEV and TEST remain unpopulated, so this package is not a fit-ready dataset and cannot support model selection or a generalization claim.
