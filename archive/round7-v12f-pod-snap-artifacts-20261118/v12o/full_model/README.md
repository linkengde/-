# v12o full-model validation outputs

All screens use the shared 458,248-atom, 200 x 200 x 300 Å model archived under `../v12n/full_model/Ag_Ti3SiC2_bicontinuous_200x200x300.data.gz`. Force snapshots and cold endpoints are gzip-compressed here. Inputs, logs, stdout, metadata, and JSON audits are retained alongside.

POD and SNAP both completed their 10-step (10 fs) no-laser cold screens without atom loss. Neither is cold-stable: see each `cold_screen_summary.json`. These results do not authorize heating.
