# F DFT result archive layout

Create this layout only after an approved structure, task-specific release, and durable writable `DFT_WORK_ROOT` are available. The current workspace overlay and absent DFT root are not an archive destination.

```text
${DFT_WORK_ROOT}/f/<task_id>/
├── input/
│   ├── structure.extxyz             # exact reviewed input; never rewrite in-place
│   ├── input_manifest.json          # source lineage, geometry hash, frozen method/release
│   └── source_manifest.json         # original source and transformation hashes
├── checkpoint/
│   └── state.gpw                    # live/recovery state; private to this F task
├── logs/
│   ├── gpaw.log
│   ├── launcher.log
│   └── progress.json
├── output/
│   ├── labeled.extxyz               # exact input geometry plus REF_forces and label metadata
│   ├── summary.json                 # RESULT_FORMAT.template.json with measured values
│   ├── verification.json            # machine-checkable checks and failures
│   └── SHA256SUMS.txt               # hashes of data files, excluding this file and archive_manifest.json
└── archive_manifest.json            # final task index; hashes data files and SHA256SUMS.txt
```

## Archive rules

- Keep input, result, and report hashes distinct. Hash raw input bytes and the ordered symbols, atom IDs, cell and PBC representations used by the verifier.
- Keep both GPAW native energy and force-consistent free energy. Name the selected target key explicitly. For finite-smearing thermal candidates, follow C's frozen label convention; do not place the two energy definitions in one unlabeled field.
- Store all per-atom force vectors in `labeled.extxyz` as `REF_forces` in eV/Å. Preserve the original symbols/order/IDs/coordinates/cell/PBC and verify them exactly against input.
- Store the raw symmetric 3×3 stress in eV/Å³ and Voigt `[xx, yy, zz, yz, xz, xy]`. For slabs, mark the vacuum-diluted 3D tensor as a raw cell stress; do not report it as bulk stress.
- Include SCF convergence, iteration count, code/runtime versions, PAW hashes, functional, cutoff, smearing, k mesh, thresholds, MPI ranks, wall time, peak memory, data role, release SHA, and verifier result.
- Keep checkpoints in the task-private recovery area. Do not push an active `.gpw` checkpoint or expose it as a finished result. If a checkpoint is needed across sessions, verify the destination and restore test before relying on it.
- The final archive is immutable after verification. `SHA256SUMS.txt` lists every payload file except itself and `archive_manifest.json`; the archive manifest records those payload hashes and the SHA-256 of `SHA256SUMS.txt`, but does not hash itself. Rerun the archive verifier from a clean read and record the verifier's own SHA/version. Any correction creates a new archive revision with an explicit parent hash.
- No DFT output archive has been created by this preparation.
