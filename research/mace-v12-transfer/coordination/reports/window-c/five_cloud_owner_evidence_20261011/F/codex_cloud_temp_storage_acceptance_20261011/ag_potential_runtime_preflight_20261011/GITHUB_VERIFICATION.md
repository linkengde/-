# F branch upload and remote verification

Repository: `linkengde/-` (public)  
Branch: `cloud-f-b1-dft`

- Report upload commit: `7f0d65cc02144ae3ae081ae5415027934a569269` (normal fast-forward push; prior head was `4addb5c3ced736589737a005d5c578331291c8d4`).
- A fresh temporary Git repository fetched that branch and commit from GitHub.
- All six report payload files were downloaded from the fetched commit; `sha256sum -c SHA256SUMS.txt` passed for each.
- The previously qualified non-science `RECOVERY_TEST.bin` fixture was fetched again from the same remote branch. Its 16,415-byte restored copy matched SHA-256 `24f3c1e665f4e0f74ea736d493704b7756871e238d2956761b001aff62921c6b`.
- The branch diff contains no `Ag.eam.alloy`, PLT parameter, or PDF payload. The local potential SHA was independently checked after remote verification.

Result: branch write, fresh remote fetch, report checksum verification, existing fixture restoration, and potential-file exclusion all PASS. The preflight itself remains HOLD solely because the signed C config does not release F.
