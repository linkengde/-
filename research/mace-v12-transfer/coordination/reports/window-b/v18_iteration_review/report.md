# V9–V18 cumulative review: clean30 versus provisional43

## Matched same-development results
|Interface|Vector V17 → V18|Separation V17 → V18|
|---|---|---|
|AgC|0.128254 → 0.135904|0.094796 → 0.213996|
|AgSi|0.061294 → 0.057137|0.190926 → 0.190539|
|AgTi|0.094163 → 0.095893|0.213471 → 0.228859|

AgC vector/max/separation worsen and loses its previous separation pass. AgSi vector/max and energy improve, but separation remains near0.19 and fails. AgTi aggregate/separation worsen slightly while maximum improves slightly. All vector gates fail; all separation gates fail; energies still pass. Clean provenance is beneficial bookkeeping, not automatic force accuracy.

### AgC
V17: shell=0.152237, outside=0.054725, longitudinal=0.068021, transverse=0.108730, pairAg=0.108179, pairX=0.202976, pairdifference=0.094796 eV/A.
V18: shell=0.161996, outside=0.054077, longitudinal=0.077872, transverse=0.111381, pairAg=0.082543, pairX=0.296539, pairdifference=0.213996 eV/A.
### AgSi
V17: shell=0.067562, outside=0.050051, longitudinal=0.040250, transverse=0.046227, pairAg=0.077038, pairX=-0.113888, pairdifference=-0.190926 eV/A.
V18: shell=0.063381, outside=0.045806, longitudinal=0.037997, transverse=0.042672, pairAg=0.086595, pairX=-0.103944, pairdifference=-0.190539 eV/A.
### AgTi
V17: shell=0.115855, outside=0.060964, longitudinal=0.062256, transverse=0.070646, pairAg=-0.018657, pairX=-0.232128, pairdifference=-0.213471 eV/A.
V18: shell=0.119216, outside=0.059365, longitudinal=0.064578, transverse=0.070888, pairAg=-0.010684, pairX=-0.239543, pairdifference=-0.228859 eV/A.

## Completion and controlled comparison limits
Model/data/reference/selection hashes match. Independently checked all99 exported state tensors equal the frozen epoch79 checkpoint; no inference was run. Complete80-epoch/Done log exists, actual threads4. Original launcher exit receipt is unavailable: training_exit_code null is preserved, not replaced with inferred zero. Numerical epoch79 artifact identity does not recover historical process exit status.
Prior V9–V17 reviews remain valid with their comparable-set limitations. V18 changes43->30 membership by removing14 partial-source rows and adding exactGamma frame0. Coverage removal, reference distribution and mesh change are inseparable in this comparison; do not identify provenance as a sole cause. Defaults retain medium/seed45/batch1/lr1e-4/wd5e-7/80epochs/energy100/force1000 and fixed final79. Reused development references are not fresh validation.

## Minimum evidence-driven V19 follow-up
A may review the six already-designed signed residual candidates, then authorize only those nonredundant acquisitions needed around AgC transverseTi and AgTi markedTi/frameworkC differential response, with AgSi pair differential as secondary coverage. Treat reused-development parents explicitly as acquisitions and freeze a fresh independent validation family before fitting. Hold recipe/checkpoint policy fixed and change only approved training labels in the next comparison; no weight/epoch changes justified by these observations alone. If no new label is authorized, report blocked coverage rather than lower gates.
Keep all current gates and sealed test labels outside training/choice. Contact screening and parent near-correlation require A review before assigning DFT owners. This recommendation is not an execution authorization, and no MD/TTM is supported.
