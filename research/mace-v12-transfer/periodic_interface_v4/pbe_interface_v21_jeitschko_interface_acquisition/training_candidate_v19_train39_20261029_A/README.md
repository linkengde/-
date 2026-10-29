# V19 train39 Jeitschko interface experiment

This is an experiment-only refit built from the frozen V19 `train36` plus three independently sourced Ag/Ti3SiC2 interface labels. It does not overwrite the approved train36/dev3 files or the previous V19 model. The three added structures share the Jeitschko–Nowotny 1967 parent family; the frozen dev3 remains the separate COD 9009647/Kisi 1998 family.

The data package contains 39 training frames and an exact byte-copy of the frozen 3-frame dev set. The original train36 bytes are the exact prefix of the train39 file. For each added frame, native extrapolated energy is mapped to `REF_energy`, GPAW forces at sigma 0.10 eV are mapped to `REF_forces`, and the free energy remains separately recorded. The added labels match dev3's slab settings and differ from train36's Gamma/fully-periodic settings.

The composition matrix remains rank 4. Its raw singular-value condition number improves from 384.76 for train36 to 40.71 after adding the three structures, because the new labels add the Ag4Si8Ti24C16 composition represented in dev3.

The experiment uses the fixed V19 recipe: seed 45, CPU, four threads, batch size 1, 80 epochs, learning rate 1e-4, weight decay 5e-7, energy weight 100, force weight 1000, with the same foundation model and `REF_energy`/`REF_forces` keys. It evaluates only the frozen dev3; no test or sealed data are included.

Status: training is running. This candidate is not an approved official dataset, and the old V19 model remains unchanged. The results must pass epoch/checkpoint verification and the existing development force and separation-force gates before any next decision.
