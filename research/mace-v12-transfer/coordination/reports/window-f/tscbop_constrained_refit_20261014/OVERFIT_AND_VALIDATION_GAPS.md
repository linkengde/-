# Overfit risk and independent-validation gaps

- Four displaced records are signed perturbations of one parent bulk geometry; the three bulk mesh records are that same geometry. Counting them as independent structures would overweight one family.
- Only one Ti site and one C site are displaced, only along z. There are no independent sites or x/y responses in the audited labels.
- There are no stress labels or converged volume/shear response records, so elastic and lattice-response behavior cannot be constrained.
- The E report's 114-atom surface diagnostic is interface-derived and frozen. Its raw input and label files are absent from this checkout, so it cannot currently be used for a direct fit or held-out score.
- No source-independent validation structure was found among the authorized, locally inspected data. The same geometry at different k meshes is not validation.
- No temperature, defect, clean-surface, or separately sourced bulk environments were assessed in this task.
- No parameter sensitivity or fit was executed. Therefore no candidate overfit test exists; any future in-sample reduction alone would be insufficient evidence.

`INDEPENDENT_VALIDATION_PASS=NOT_TESTED`. A later validation set should contain at least an independently sourced bulk cell or distinct displacement family, plus a separately frozen strain/surface family appropriate to the solid-skeleton use. Its labels must remain unopened during fitting.
