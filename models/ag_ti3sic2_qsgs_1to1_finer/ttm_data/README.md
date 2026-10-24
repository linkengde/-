# Ag TTM candidate table

`Ag_TTM_hybrid_figure_fit_300_20000K.csv` is copied byte-for-byte from the user-provided `Ag_TTM_processed_20261024.zip` member with the same name.

- SHA-256: `9cb40987150ce0029411b64671b5dd5999ae2f0b6ea7e173fe214f1130fb8a5b`
- Selected columns: `Te_K`, `Ce_eV_A-3_K-1`, `G_eV_ps-1_A-3_K-1`, `ke_eV_ps-1_A-1_K-1`, and `Ue_delta_from_300K_eV_A-3`.
- Range: 300–20,000 K at 100 K spacing.
- These are cross-paper digitized/fitted candidates (Lin et al. for `Ce`/`G`, Karna & Giri for `ke`), not a production-validated self-consistent Ag parameterization.
- The reference solver must stop outside this temperature range; it does not extrapolate or clamp.
- An optional low-temperature closure in the reference solver extends 1–300 K with `Ce(T)=Ce(300)·T/300`, `G(T)=G(300)`, and `ke(T)=ke(300)`, and integrates this `Ce` to obtain `Ue`. This is an explicit, unvalidated Sommerfeld/constant-coefficient assumption, not part of the supplied table. It is opt-in and must be sensitivity-tested before use.

The source archive’s README and audit JSON remain in the user attachment. This copy is only for reproducible, read-only parameter access by the reference code.
