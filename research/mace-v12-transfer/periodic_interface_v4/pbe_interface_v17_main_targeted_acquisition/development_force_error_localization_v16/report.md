# V16 force-error localization on V17 development structures

This is a diagnostic of the frozen V16 model on the three authorized V17 development structures. The labels remain validation data for V17; this report does not move them into training.

The local shell is atoms within 4.5 Å MIC of either marked Ag–X atom. Shell RMSE is the root mean square of per-atom force-error vector norms. Directional components are projected onto the marked Ag-to-X axis. These diagnostics identify where new data may help; they do not establish the cause by themselves.

| Interface | Pair distance (Å) | Vector RMSE | Separation error | 4.5 Å shell RMSE | Outside-shell RMSE | Largest residual |
|---|---:|---:|---:|---:|---:|---|
| AgC | 2.660 | 0.1376 | 0.2044 | 0.1632 | 0.0596 | Ti#604600 (0.4590) |
| AgSi | 2.888 | 0.0729 | 0.2035 | 0.0805 | 0.0591 | Si#366936 (0.2035) |
| AgTi | 3.223 | 0.1015 | 0.1856 | 0.1253 | 0.0650 | Ti#479027 (0.3154) |

## What the residuals point to

- **Ag–C:** the largest error is Ti#604600 (0.4590 eV/Å), with a 0.4542 eV/Å component transverse to the marked Ag–C axis. Its local neighbors include Ag#751743 at 2.154 Å and several framework C atoms at about 2.10–2.22 Å. The marked C#604804 also has 0.3024 eV/Å error. The local shell RMSE is much larger than outside the shell.
- **Ag–Si:** the largest error is Si#366936 (0.2035 eV/Å); marked Si#366934 and nearby Ag#685172 also have predominantly transverse residuals. Local shell error is only moderately above the outside shell, so this is not confined to the marked pair alone.
- **Ag–Ti:** marked Ti#479027 has 0.3154 eV/Å error and nearby framework C#480805 has 0.2829 eV/Å. The Ti–C distance between them is 1.869 Å. The local shell RMSE is roughly twice the outside-shell value.

## Data decision

Prioritize the already verified eight V17 positive/negative training diagnostics: their Ag–C, Ag–Si transverse/framework, and Ag–Ti framework moves target the same residual families. Keep these three development structures out of training for the first V17 comparison. If V17 still misses the gates, add paired DFT perturbations of the implicated framework-neighbor environments and transverse registry, then reserve a new untouched family for validation. Do not tune force weight before checking this data-coverage change with the V16 recipe fixed.

The six withheld V17 structures and their labels were not accessed.
