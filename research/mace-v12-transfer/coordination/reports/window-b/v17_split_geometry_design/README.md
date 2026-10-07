# V17 split geometry design — Window B

The final candidate bundle is in [final_set/](final_set/README.md). It contains three development-validation and six withheld-test geometries, with source provenance, per-ID displacement records, duplicate/contact screening, and SHA256 manifests.

During screening, the three previously unassigned V16 geometry proposals were found to be exact duplicates of the frozen V16 blind-holdout inputs. They were excluded. Geometry-only evidence is preserved under [attempt_diagnostics/attempt_01/](attempt_diagnostics/attempt_01/diagnostic.json); no holdout labels or calculation outputs were read.

The final set has no exact/near hits against the 197 historical model-split frames or 38 published DFT input geometries, and its nine role geometries are pairwise separated under the documented criteria. These are geometry-triage results only. All candidates inherit current small-cluster parent motifs; A must review the roles before training, and withheld-test labels must stay outside checkpoint/parameter selection.

See [final_set/audit_report.md](final_set/audit_report.md) for the full report, [final_set/role_manifest.json](final_set/role_manifest.json) for assignments and hashes, and [final_set/screening_manifest.json](final_set/screening_manifest.json) for machine-readable screening.
