# V16 controlled data-augmentation screening

35 training frames = V15 retained29 + six verified distance/registry acquisitions.
Training settings unchanged: foundation, seed45, CPU, batch1,80 epochs, lr1e-4,
energy100/force1000. Historical valid2/test3 are byte-preserved regression data.
Three new frozen local-registry probes excluded from training/validation/model
selection. Their small-cluster motif correlation and inherited reference-method
uncertainty remain; passing these gates cannot establish broad model validity.
Old V15 scored probes are history, not fresh V16 validation.

Build with existing GPAW Python; preflight and train use reviewed MACE environment.
run_v16_training.sh requires explicit repo/output directory; refuses existing output
and active MPI. Freeze model selection/completion evidence before evaluate_v16.py.
No long MD/TTM authorization follows numerical training or narrow screen results.
