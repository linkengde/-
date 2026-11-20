# v12p training candidate

Parent v12o; retains the 128 train / 12 validation split and Si01 center 104359 as a development holdout. Adds four converged DFT force-label crops from the v12o SNAP 10 fs no-laser collapse: two Ag-Si, one Ag-C, one Ag-Ti. Only central atom pair force rows use DFT; all crop-edge forces use the POD v10 anchor. DFT labels are not equilibrium states and do not by themselves establish transferability. Fit POD/SNAP, check holdout and full-model static/cold screens before any heating decision.
