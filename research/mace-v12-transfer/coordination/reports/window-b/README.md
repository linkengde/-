# Window B analysis reports

Window B completed its four assigned DFT labels, the v13 per-atom residual report, and one v14 training/evaluation cycle. V14 screened `SCREEN_FAIL`; no production MD/TTM is authorized.

## Current assignment

While A runs three residual-targeted v15 candidate DFT labels, B analyzes atom-level force errors on the three scored v14 blind holdouts. The report identifies residual atoms, force directions, contact-local and remote-region errors, and nearest neighbors. It is post-hoc diagnosis only; it must not be used to change the completed v14 model or its data.

Start with `coordination/START_WINDOW_B.md` and run `run_v14_force_localization.sh` after pulling `main`. The script verifies the exact checkpoint and DFT reference hashes, writes `v14_force_localization.json` and `.md`, and publishes B task progress.

Window A owns the v15 DFT queue, main work log and global manifest. B must wait for A's next explicit assignment before another training cycle. Any next model needs fresh blind holdouts.
