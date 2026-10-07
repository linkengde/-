# V17 provisional run03 result

The corrected CPU run completed 80 epochs and exported the predeclared final epoch 79. The runner reached its success receipt; its log records all epochs 0–79, `Training complete`, `Done`, and `Loaded Stage one model from epoch 79 for evaluation`. Four CPU threads were used. No test file or withheld labels were opened.

| Development structure | Energy abs. error (meV/atom) | Force-vector RMSE (eV/Å) | Separation-force abs. error (eV/Å) | Result |
|---|---:|---:|---:|---|
| AgC | 2.755 | 0.1283 | 0.0948 | Energy and separation pass; force fails |
| AgSi | 1.771 | 0.0613 | 0.1909 | Energy passes; force and separation fail |
| AgTi | 1.221 | 0.0942 | 0.2135 | Energy passes; force and separation fail |

Provisional gates are 10 meV/atom energy, 0.05 eV/Å force-vector RMSE and 0.10 eV/Å separating-force error. Overall development status is **FAIL**. The three structures were already used to diagnose V16; this is not independent blind validation. All eight signed training diagnostics and the inherited data-source limitations remain relevant. Do not use this model for production MD/TTM.

Full per-structure evidence is in `../evaluation_dev_01/`; model, final checkpoint, complete logs, training metrics, selected data and scripts are enumerated in `SHA256SUMS.txt`.
