# V17 energy convention and parent geometry hash audit

## Scope

This audit traces three V14 **training** rows and their source GPAW archives, then
checks the parent relation and reference convention for the eight published V17
paired training diagnostics. It reads no V17 development-validation output and no
withheld-test label, summary, log, or score. No DFT, MACE inference/training,
production-data edit, MD, or TTM was performed.

The V14 source dataset is
`periodic_interface_v4/mace_periodic_v14_interface_energy/data/train.extxyz`,
SHA256 `3323cdbd3081fc8d91673662c5f83c8da93ba3b4d7a23428f6b70d0e7a4c9faf`.
All archive input/output/summary/log hashes and the eight V17 diagnostic archive
hashes were rechecked. The complete machine evidence is in `audit.json`; the
compact per-parent table is `parent_evidence.csv`.

## Parent hash trace

The three upstream `source_frame_geometry_sha256` strings first appear in the
committed `v15_registry_candidates/candidate_screening_manifest.json` added by
commit `a3b1299`. That commit contains the three candidate extxyz files, an audit
report, and the large manifest, but no candidate-generation or hashing script.
The manifest has no hash-schema declaration. The V14 acquisition manifest records
input-file hashes, and its builder copies energies and forces; neither defines
these frame-geometry hashes. Later V16 proposal code copies the field from the
manifest, and the V17 split generator copies `source_parent.frame_geometry_sha256`
without recalculating it. The original serialization therefore cannot be
reproduced from committed evidence.

The reproducible V17 canonical function is documented in
`v17_split_geometry_design/final_set/generate_split_candidates.py::geometry_sha`
and independently repeated in `audit.py`. It hashes UTF-8 bytes of compact,
key-sorted JSON (`sort_keys=True`, `separators=(',', ':')`) with these fields:

- `symbols`, in extxyz row order;
- `positions_A`, Cartesian ASE coordinates in Å rounded to 10 decimals;
- `cell_A`, the ASE cell in Å rounded to 10 decimals;
- boolean `pbc` values;
- `lammps_id`, in the same row order.

SHA256 of those bytes gives:

| Parent | V14 frame | Declared legacy hash | V17 canonical hash | Exact match |
|---|---:|---|---|---|
| AgC `AgC_registry_strain_acq_v14_01` | 29 | `72d649542158aade39a3cf50dc73759ad25df8a7ff7b2b117f1aa0c6ffc513a3` | `76fdf4d13972b2ebed7b0476f5166607f003ffa690b36d53707afec196aaccf8` | No |
| AgSi `AgSi_registry_strain_acq_v14_01` | 30 | `4a5540c7fd824a01df8fcb6155e62c7ff045ec36d637ac5f405247c9b5dbe4e5` | `15d350743e243e83f26abb27101e983b8970bb63a1569d56b6b65d56fcd43226` | No |
| AgTi `AgTi_registry_probe_v13_01` | 28 | `18dc1be733a506b748386175293d2f094e720e3e1c9dc96e259623088c005b7c` | `0ad9bc5d0d6636b02c31048d0c33aa198786e5d4fa784a668701c1e744bcd764` | No |

The hash mismatch is a metadata/provenance gap. The three dataset rows still match
the specified V14 dataset SHA, frame index, config type, composition, atom count,
and persistent IDs. Their archived DFT outputs match row order, elements,
coordinates, cell, and PBC. For all eight V17 diagnostic inputs, the recorded
source row, atom ID/order and element mapping, cell, and PBC match; the moved IDs
and signed displacements match the input manifest. Maximum changed-atom vector
discrepancy is `4.24e-9 Å`; maximum undeclared displacement is `0 Å`. All eight
input/output archives independently verify. Thus the old digest mismatch does not
invalidate the DFT archives or their parent association, but the old hash itself
must remain marked as unreproduced.

## V14 energy and force lineage

The V14 runner constructs GPAW with PBE, plane waves at 500 eV, Gamma-only
`(1,1,1)` k-points, Fermi-Dirac width 0.1 eV, mixer beta 0.05 / 8 old densities /
weight 100, and `1e-5` energy, density, and eigenstate convergence thresholds.
It calls `atoms.get_potential_energy()` with the default `force_consistent=False`
and `atoms.get_forces()`. The archive summaries record GPAW 26.7.0 and four MPI
ranks; the GPAW logs contain the Fermi-Dirac width and completed SCF step counts.
The V14 dataset builder copies `PW_PBE_energy_eV` to `REF_energy` and
`PW_PBE_forces` to `REF_forces` without transformation.

| Interface | Frame / atoms | `REF_energy` = archive native energy (eV/cell) | Log `Extrapolated` (eV, rounded) | Log `Free energy` (eV, rounded) | SCF steps / ranks |
|---|---:|---:|---:|---:|---:|
| AgC | 29 / 48 | -257.5833982481984 | -257.583398 | -258.740849 | 58 / 4 |
| AgSi | 30 / 51 | -208.69345818474375 | -208.693458 | -209.631572 | 49 / 4 |
| AgTi | 28 / 49 | -176.57984627488767 | -176.579846 | -177.224883 | 38 / 4 |

For each row, `REF_energy`, the source extxyz `PW_PBE_energy_eV`, the source
summary `energy_eV_cell`, and the V14 dataset-manifest energy agree exactly. They
match the GPAW log's `Extrapolated` line to its six-decimal precision and do not
match its separate `Free energy` line. The log-only free energies are rounded to
six decimals; the V14 output and summary did not preserve their full precision.
The V14 frame's `REF_forces` is an exact elementwise array copy of the source
archive's `PW_PBE_forces`. It is the direct `get_forces()` result for the logged
0.1 eV Fermi-Dirac calculation.

The historical `source_method` string says “GPAW 26.7.0 PW-PBE 500 eV Gamma
single point”; it omits smearing, the native/free energy distinction, and an
explicit force/energy convention record. The archived script, summary, and logs
fill in the method and energy lineage. They do not make the legacy row metadata
complete by themselves.

## V17 paired diagnostic and later provenance evidence

All eight paired diagnostic archives report `PASS`; the audit rechecked their
input and archive hashes, four ranks, convergence, finite labels, structure
identity, and method. Their summaries retain `energy_eV_cell` and
`free_energy_eV_cell` separately and explicitly state “GPAW native extrapolated
energy, force_consistent=False; free_energy recorded separately.” The method is
the same PBE/PW-500/Gamma/Fermi-Dirac-0.1/GPAW-26.7 setup as the three V14 source
archives.

The already-published paired-response report found that the free-energy central
secant had a smaller residual against the projected DFT force response for three
of four finite-displacement pairs (AgSi transverse, AgSi framework neighbor,
AgTi framework neighbor); AgC was close under both definitions. These secants are
finite-displacement diagnostics, not proof that an archive or force label is
incorrect. The evidence supports preserving both energy fields and recording the
chosen training target explicitly.

The V16 dataset-integration audit gives a useful later schema example: all six
new V16 archive additions explicitly declare native extrapolated energy with
`force_consistent=False` and retain free energy separately. That audit still
marks inherited physical-method consistency as `UNKNOWN`; its six new records do
not retroactively resolve older rows or the legacy geometry-hash algorithm.

## Recommendation for V17 integration

1. **No relabeling is required for these three parent energy values** if V17 keeps
   the established native-extrapolated cell-energy target. Use the existing
   `REF_energy` for these parents; use each new V17 summary's `energy_eV_cell` for
   new `REF_energy`, while retaining `free_energy_eV_cell` in a separate field.
2. Keep the original force arrays unchanged. Record the force extraction,
   smearing, code version, cutoff, k-points, convergence settings, archive hash,
   and the energy convention alongside each frame. Do not state that the old
   `REF_forces` are force-consistent with the extrapolated energy: that flag was
   not persisted in the legacy frame, and the paired secants show a finite-smearing
   convention difference worth tracking.
3. Preserve each declared parent hash as a legacy provenance value. Add the
   reproducible V17 canonical hash in a separate field and mark the legacy hash
   algorithm unknown. Do not overwrite source metadata or describe this mismatch
   as an archive failure.
4. This audit resolves only the three assigned V14 parents. The V16 audit still
   reports method consistency as unknown across inherited training records; a
   corpus-wide claim of uniform DFT method/energy provenance needs the separate
   inherited-label audit. If A chooses free energy as the training target instead,
   do not silently substitute the six-decimal V14 log values; first set a
   precision-preserving conversion or relabeling policy for every included row.

## Reproduction files

`audit.py` regenerated `audit.json` and `parent_evidence.csv` with the command in
`README.md`. `SHA256SUMS.txt` covers all published deliverables except itself.
