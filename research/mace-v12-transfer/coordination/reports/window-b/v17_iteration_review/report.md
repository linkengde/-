# V9–V17 cumulative review and matched V16/V17 development comparison

|Interface|Vector RMSE V16 → V17|Separation error V16 → V17|
|---|---|---|
|AgC|0.13764 → 0.12825|0.20439 → 0.09480|
|AgSi|0.07289 → 0.06129|0.20347 → 0.19093|
|AgTi|0.10153 → 0.09416|0.18565 → 0.21347|

All three vector RMSEs improve but remain above0.05 eV/A. AgC separation improves and passes0.10; AgSi improves but fails; AgTi separation worsens and fails. AgTi maximum atom error also worsens despite lower aggregate RMSE. All energies improve and pass10 meV/atom. Overall remains FAIL. References/model/data/epoch79 selection hashes checked exactly. These are reused development structures, not fresh blind tests.

## Cumulative evidence and limits
Reuse the hash-pinned V9–V16 report: local transverse/framework coverage remained the leading actionable hypothesis, provenance uncertainty a separate risk, historical roles/metrics limit comparisons, and no controlled optimizer ablation proved a parameter cause. V17 adds eight signed local diagnostics; matched aggregate improvements support a useful coverage contribution but do not isolate causation. Selection policy changed to fixed epoch79 and earlier runs had checkpoint-export issues, so this is not a perfectly isolated data-only experiment. The earlier claim of12 additional declaration-only rows is superseded by original archive reconciliation: current pool has29 direct and14 partial.

V16 localization identifies AgC Ti604600 transverse response, AgSi framework/Ag transverse neighbors, and AgTi Ti479027–C480805 shell response. V17 published evaluation contains only aggregate/max metrics, no per-atom prediction vectors. Therefore V17 spatial, transverse and framework changes cannot be reconstructed without forbidden new inference; do not transplant V16 atom identities onto V17 maxima.

## V18 recommendations, ordered by evidence
1. Maintain gates and record AgTi separation/max regression as an explicit target. Obtain/publish V17 development per-atom diagnostics from existing A artifacts if available before choosing perturbation directions. Do not infer one cause from aggregate errors.
2. Prepare the directly supported30-frame proposal (core29 plus exact frame0 relabel) and let A decide a provenance/data-membership comparison with fixed recipe and the same development references. This removes14 partial frames but changes coverage and adds a mesh-changed relabel, so an error difference cannot establish label inconsistency alone. Gamma convergence of short frame0 directions is unresolved.
3. If coverage-driven augmentation is chosen, prioritize paired AgTi marked Ti/framework-C and AgC transverse Ti-shell environments; keep prior paired signs and new training roles explicit. Require fresh held-out family after development tuning, without opening current sealed tests for choice.
4. No evidence presently identifies a better learning rate or force weight. Preserve seed45, medium foundation, batch1, lr1e-4, wd5e-7, energy100/force1000 and80 epochs for the minimum next controlled comparison. Do not vary data, weights and epoch policy together.

Minimum controlled follow-up for A review: one clean30-versus-provisional43 fixed-recipe comparison, reporting vector/max/separation on the same development set, provenance/mesh caveats and coverage removals; not an authorized training launch. No claim of thermal/morphology applicability and no MD/TTM.
