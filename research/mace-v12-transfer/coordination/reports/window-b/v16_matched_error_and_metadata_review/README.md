# V16 matched error and metadata review

This is a retrospective V15/V16 comparison on the same three previously scored V15 holdouts and the same three previously scored V16 holdouts. The selected checkpoints are frozen; the results are not used for model selection or as fresh blind validation.

Run from the repository root with the existing MACE CPU environment:

```bash
MPLCONFIGDIR=/tmp/mpl-config XDG_CACHE_HOME=/tmp/xdg-cache \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/workspace/.venvs/mace-v14-assist/bin/python -B \
research/mace-v12-transfer/coordination/reports/window-b/v16_matched_error_and_metadata_review/analyze_matched_v15_v16.py \
--repo-root /workspace/- \
--output-dir /workspace/-/research/mace-v12-transfer/coordination/reports/window-b/v16_matched_error_and_metadata_review/reproduction_01
```

The script verifies V16 metadata pins and unchanged data splits, frozen model and evaluation hashes, archived reference/input hash chains and score reproduction before saving predictions. It refuses a nonempty output directory. The verified run is `analysis_04/`; `analysis_01/` through `analysis_03/` retain failed diagnostic runs from fixing the analysis script and did not produce model-selection results.

`SHA256SUMS.txt` covers all files in this report directory except itself.
