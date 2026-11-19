# v12n POD/SNAP candidate

Parent: v12m. This candidate preserves the existing 128 training / 12 validation split and held-out Si01 center 104359. It adds two GPAW PBE/PW350 Gamma DFT force labels for the same C-C pair at 1.386 and 1.114 Å during a cold-run collapse precursor. Only the central two C force rows per crop are DFT targets; edge rows use v10 POD force anchors.

This is a training candidate only. It must pass the independent Si01 holdout, full-model static force checks, and full-model cold stability before any heating trial.
