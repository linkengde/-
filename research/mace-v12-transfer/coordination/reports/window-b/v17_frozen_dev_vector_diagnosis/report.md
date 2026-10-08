# Frozen V17 same-development vector diagnosis

Published V17 energy/vector/separation metrics reproduced within1e-4 meV/atom energy and1e-6 eV/A force (CPU/runtime arithmetic tolerance). Only approved frozen epoch79 model and three scored development references were used; no unchanged V16 inference, fitting, DFT or withheld access. Signed separation error is projected X error minus projected Ag error, so improvements in bulk RMS need not improve this two-atom difference.

## AgC_v17_dev_validation_01
V16: vector=0.137644, shell=0.163219, outside=0.059640, transverse=0.114889; Ag signed=0.087898, X signed=0.292288, separation signed=0.204390 eV/A.
V17: vector=0.128254, shell=0.152237, outside=0.054725, transverse=0.108730; Ag signed=0.108179, X signed=0.202976, separation signed=0.094796 eV/A.
Largest V17 residual IDs: [(604600, 'Ti', 0.42846), (751743, 'Ag', 0.30008), (604805, 'Ti', 0.27738), (688433, 'Ag', 0.27473), (604603, 'C', 0.24528)]
## AgSi_v17_dev_validation_01
V16: vector=0.072894, shell=0.080541, outside=0.059119, transverse=0.054142; Ag signed=0.074252, X signed=-0.129221, separation signed=-0.203473 eV/A.
V17: vector=0.061294, shell=0.067562, outside=0.050051, transverse=0.046227; Ag signed=0.077038, X signed=-0.113888, separation signed=-0.190926 eV/A.
Largest V17 residual IDs: [(366934, 'Si', 0.16153), (366936, 'Si', 0.14446), (685172, 'Ag', 0.11909), (683609, 'Ag', 0.08996), (368827, 'Ti', 0.08862)]
## AgTi_v17_dev_validation_01
V16: vector=0.101534, shell=0.125250, outside=0.065030, transverse=0.072881; Ag signed=-0.017617, X signed=-0.203264, separation signed=-0.185647 eV/A.
V17: vector=0.094163, shell=0.115855, outside=0.060964, transverse=0.070646; Ag signed=-0.018657, X signed=-0.232128, separation signed=-0.213471 eV/A.
Largest V17 residual IDs: [(479027, 'Ti', 0.3482), (480805, 'C', 0.26939), (479028, 'Ti', 0.18754), (479021, 'Ti', 0.15173), (480980, 'Ti', 0.13633)]

## Acquisition recommendations for A review
First rank contact-shell neighbor residuals from this report and preserve pair contributions separately. For AgTi, target marked Ti plus nearest framework-C response and transverse Ti/Ag registry; propose opposite-sign small displacements along the observed residual directions, with geometry/contact checks and new-role isolation before any DFT approval. AgC prioritize large Ti transverse shell residual; AgSi prioritize marked Si/Ag differential projection and neighboring framework response. These directions are diagnostic hypotheses, not evidence of a sole cause or authorization to calculate.
The existing eight pairs probe particular parent environments; residual maxima on moved development registries can require additional neighbor coverage. A fixed-recipe clean30 comparison is useful but combines provenance filtering, coverage removal and mesh-changed frame0; separate these confounds. Keep force/separation gates unchanged and sealed tests unavailable for tuning.
