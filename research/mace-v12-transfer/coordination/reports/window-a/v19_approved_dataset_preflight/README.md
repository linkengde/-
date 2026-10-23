# V19 approved dataset and preflight

A reviewed B audit passed with a shared-source limitation. Window A approved and assembled `A_APPROVED_V19_DATASET` with train36/dev3/test0. The train36 file is byte-identical to the prior candidate; `valid.extxyz` contains the archive-verified AgC, AgSi, and AgTi development labels.

The read-only B audit script was rerun, its report checksums passed, and `preflight_v19.py` returned `V19_PREFLIGHT_PASS`. The locked runtime is Python 3.12.14, CPU PyTorch 2.5.1, mace-torch 0.3.16, and ASE 3.29.0; `pip check` reported no broken requirements.

The three dev structures are distinct under the audited geometry screen but share one Ti3SiC2 COD 9009647 parent source family; the upstream source bytes were not externally verified. Slab k-point, smearing, and vacuum convergence remain open. The dataset has no test file, and no sealed/test labels were opened.

`metrics.json` records the data and code hashes, environment lock, fixed 80-epoch CPU recipe, preflight result, and disk availability before training. Training has not yet started at the time of this status record.
