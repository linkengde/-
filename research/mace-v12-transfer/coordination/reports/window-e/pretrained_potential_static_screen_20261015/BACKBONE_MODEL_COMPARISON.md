# Ti₃SiC₂ backbone comparison

All comparisons use the same archived DFT coordinates and atom mapping. Force metric is `sqrt(mean_i(||F_model_i-F_DFT_i||²))` in eV/Å. These are screening references, not blind validation.

| Model | 48-atom bulk RMSE | 114-atom TSC-only RMSE | TSC-only Ti | Si | C | Assessment |
|---|---:|---:|---:|---:|---:|---|
| MACE-MP-0b3-medium | 0.03817 | 0.50399 | 0.53707 | 0.77422 | 0.17213 | Best 48-atom baseline; larger slab error. |
| MACE V14 | not available in cited 48-atom per-atom report | 0.35214 | 0.39829 | 0.42477 | 0.20936 | 135-derived comparison only. |
| MACE V17 | not available in cited 48-atom per-atom report | 0.30694 | 0.36215 | 0.38243 | 0.11306 | Best listed existing MACE TSC-slab RMSE; bulk not reported here. |
| MACE V18 | 2.88885 | 2.15794 | 2.18985 | 0.41129 | 2.58694 | Severe known pure-backbone regression; do not overlook. |
| SevenNet-MF-0 | 0.05676 | 0.40616 | 0.41391 | 0.63549 | 0.18973 | Worse than MACE-MP bulk baseline and V17 slab; no consistent improvement. |
| MatterSim-v1.0.0-5M | 0.04562 | 0.20586 | 0.12598 | 0.40734 | 0.14269 | Lower TSC-slab error than compared MACE models, but bulk is slightly worse than MACE-MP; interface model-card caveat applies. |
| Old EAM+TSCBOP+Morse | not available in same 48-atom DFT audit | 2.56241 | 2.62420 | 2.60715 | 2.44326 | Retained only as historical hybrid baseline; its large error is dominated by present backbone and interface contributions. |

The 48-atom structural family is a numerical control and may overlap pretrained training distributions. The 114-atom subset derives from the same 135-atom interface parent and is not an independent family. No candidate passes a broad backbone validation gate from these data alone.
