# v12q candidate

Based on v12p, retaining the 128 training / 12 validation split and the Si01 center 104359 holdout. Adds four GPAW PBE/PW350 Gamma local force-label crops selected from short Ti-C contacts in the 10 fs SNAP v12p no-laser cold endpoint. Only the central Ti-C pair forces are DFT training targets; local edge forces use the existing POD v10 anchor. These transient high-stress crops are active-learning candidates and do not by themselves establish transferability. Run the Si01 holdout, full-model static checks, and no-laser cold stability checks before considering short heating.
