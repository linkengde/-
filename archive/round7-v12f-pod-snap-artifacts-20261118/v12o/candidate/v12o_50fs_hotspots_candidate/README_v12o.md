# v12o POD/SNAP candidate

Parent: v12n. This candidate preserves the existing 128 training / 12 validation split and held-out Si01 center 104359. It adds two GPAW PBE/PW350 Gamma DFT force labels from the new full-model cold 50 fs Ti-C and C-C short-contact environments. Only the two central target-atom force rows in each crop are DFT targets; crop-edge force rows use v10 POD anchors.

This is a training candidate only. It must pass the independent Si01 holdout, full-model static force checks, and full-model cold stability before any heating trial.
