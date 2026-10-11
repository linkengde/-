# D Ag64 step 1600 continuation release

This independent C release package authorizes only task `C-D-AGL-1250K-rho097-seed01` from global step 1600 through a hard stop at step 2000. The release is bound to the D branch head and archived step-1600 checkpoint identified in the signed JSON.

The step-1300 package remains unchanged. This directory contains a separate step-1600 adapter/core copy because the earlier adapter was pinned to the step-1300 source. No potential binary, third-party paper, or checkpoint archive is included.

The configuration retains the existing 600 s wall budget and unchanged scientific parameters. If that budget is reached before step 2000, execution must stop and preserve the latest verified checkpoint; no automatic continuation is authorized.

C performed only read-only remote archive and checkpoint verification plus release/config validation. No calculator was constructed and no MD, DFT, or training was run. This commit does not dispatch work to D.
