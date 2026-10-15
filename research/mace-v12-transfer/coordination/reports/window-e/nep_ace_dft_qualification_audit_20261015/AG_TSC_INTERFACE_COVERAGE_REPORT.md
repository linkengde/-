# Ag–Ti₃SiC₂ interface coverage

## Confirmed labelled configurations

The master inventory includes 13 confirmed Ag–Ti family, 25 Ag–Si family and 16 Ag–C family geometries. These are small periodic interface/slab configurations with four elements in the slab models; their source records include acquisition, registry, distance and residual-mode cases. Counts are unique ordered geometries after exact-source deduplication, not independent parent families. Their nearest unlike-element distance ranges are:

| Family | Ag–Ti (Å) | Ag–Si (Å) | Ag–C (Å) |
|---|---:|---:|---:|
| Ag–Ti family | 2.110–2.485 | 2.963–3.419 | 2.142–2.628 |
| Ag–Si family | 2.443–5.033 | 2.001–2.967 | 1.929–5.986 |
| Ag–C family | 2.107–2.611 | 3.788–4.389 | 1.989–2.465 |

These ranges report the nearest pair in each whole cell; they do not imply every pair type was sampled at every separation. The 135-atom low-mismatch archive contributes two full four-element geometries (initial fixed and relaxed) from the same parent family and component-only references from the relaxed parent. It is not two independent interface families. Its relaxed geometry has one local registry/termination and no broad liquid-contact sampling.

## Limitations

There is useful static interface information, so interface coverage is not zero. However, source parent IDs and method details are incomplete for part of the older archive, and the split histories reuse related acquisition/slab families. There is no verified temperature-dependent liquid Ag contacting a solid Ti₃SiC₂ surface, no broad independent set of terminations/registries, and limited matched energy/force/stress data. The 77-atom three-gap structures have no DFT labels and do not increase label coverage.

Pair-isolated or slab-family calculations can help initialize local interactions but cannot replace DFT for mixed four-element environments and should not be counted as a complete liquid-solid interface training set.

`INTERFACE_DATA_READY = INSUFFICIENT_DATA`. Reuse the verified static interface labels after metadata normalization; add a small, parent-diverse set of mixed four-element interface and liquid-contact DFT configurations.
