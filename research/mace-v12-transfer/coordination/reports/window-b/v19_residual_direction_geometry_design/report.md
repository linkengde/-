# V19 label-free residual geometry proposals

18 proposals: signed amplitudes0.03/0.06/0.10A; one normalized mode per interface. AgC Ti604600 follows its observed transverse residual direction. AgTi Ti479027 and frameworkC480805 move oppositely along marked pair axis to probe differential framework response. AgSi Si366934 and Ag685172 move oppositely along axis. Two-atom relative displacement is twice nominal amplitude; all vectors/IDs/source hashes in manifest. IDs/order/cell/PBC preserved; inputs contain geometry and pair marker only.

Recommended signed pairs: ['AgC_residual_mode_m_06_v19_proposal', 'AgC_residual_mode_p_06_v19_proposal', 'AgSi_residual_mode_m_03_v19_proposal', 'AgSi_residual_mode_p_03_v19_proposal', 'AgTi_residual_mode_m_03_v19_proposal', 'AgTi_residual_mode_p_03_v19_proposal']

Prefer0.06A if both contacts pass, else0.03A, then0.10A; larger displacement improves finite-response signal but perturbs more contacts, and is reserve not an automatic compute request. Species minima use parent-specific0.05A allowance rather than all-element1.75A threshold. All passing near hits are reported: these deliberately probe a reused development family and are highly correlated with existing geometries. They are training-acquisition proposals, not blind or morphology-independent tests. Screening compares permitted historical train/dev/acquisition/prior candidates only; sealed comparison uses exact digests, so near sealed overlap remains unknown.

A must review contact/near-correlation tradeoffs and explicitly assign DFT owners before labeling. Do not label18 redundant amplitudes automatically. No inference/DFT/training/MD/TTM or official edits performed.
