# MACE v11 energy focused trial

This is an isolated fine-tune using the exact v10 train/valid/test split. The four crystalline references and AgTi PBE interface energies remain in training; the AgTi 2.70 A frame and the external AgSi/AgC structures remain excluded.

The only training-objective change from v10 is `energy_weight=100` (v10 used 1); `forces_weight=1000`, foundation model, learning rate, data, batch size, and maximum epoch count are unchanged. This increases the per-atom energy contribution to the objective while retaining force supervision.

The split files were copied byte-for-byte from v10; `data/dataset_manifest.json` records source hashes. The v11 model is a research candidate only. It must beat v9/v10 on frozen test and external holdouts before any LAMMPS smoke or MD use.
