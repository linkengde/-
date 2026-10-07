# V17 split geometry design

This bundle proposes one development-validation geometry and two withheld-test geometries per interface. All files contain geometry and identity only; no energy/force labels or DFT results are included. The V16 proposals were rejected after exact matches to V16 blind inputs were found; see the parent report. See audit_report.md for screening rules, limitations and role boundaries. A must review the role_manifest.json before training.

Candidates live in geometries/. source_inventory.json lists hashed geometry-only inputs and historical split files. screening_manifest.json contains machine-readable duplicate/contact results. SHA256SUMS.txt covers all bundle files except itself.
