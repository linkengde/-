# V16 dataset integration audit

## Scope and result

Read-only review of A's integrated `mace_periodic_v16_interface_energy` data,
the six acquisition sources, the frozen V16 holdout input geometries, and the
training/preflight entry points. No DFT, MACE training or inference, MD, or TTM
was run. The three blind calculation directories and their energy/force values
were not opened.

The data rows and source archives pass the structural checks. One definite
manifest defect remains: provenance summary counts still describe the 29-frame
parent training set although the current training split contains 35 frames.
The new six labels themselves pass source identity, archive, method, geometry,
and finite-label checks. This is a metadata/accounting defect, not evidence of
bad labels.

## Checks that passed

- V16 splits contain **35 train / 2 valid / 3 test** frames. Every energy and
  force label is present and finite; the training composition matrix has rank 4.
- All **29 retained V15 training frames** match the first 29 V16 training rows
  in atom order, positions, cell, PBC, arrays, and frame metadata.
- V16 valid and test files are byte-identical to both the V15 and V14 files.
- The five excluded-row records remain intact: four LCAO frames and the single
  unresolved AgTi frame. The inherited provenance ledger contains the expected
  29 train + 2 valid + 3 historical test records; the six new rows are separately
  listed in `V16_additions`, as the manifest's scope note states.
- All six acquisition records map one-to-one to their input manifest and
  archived PW-PBE outputs. Input/output hashes, archive checks, four MPI ranks,
  frame IDs, cell/PBC, positions, REF energy, and forces match. All six use
  GPAW 26.7.0, PW-PBE/500 eV, 1x1x1 Gamma, 0.1 eV smearing, and the same recorded
  native extrapolated-energy convention (`force_consistent=False`; free energy
  recorded separately). The V16 source ledger pins the source summaries and
  archives.
- The three frozen holdout **inputs** have the expected hashes and atom IDs.
  None exactly or nearly overlaps any of the 40 V16 train/valid/test geometries.
  The geometry screen's exact/near classification is independent of PBC; this
  was checked by toggling PBC on an in-memory copy of each input. No holdout
  calculation output or reference label was read.
- The V16 launcher retains the V12 foundation model, seed 45, CPU, batch size 1,
  80 epochs, learning rate `1e-4`, weight decay `5e-7`, energy weight 100, and
  force weight 1000. It passes only the V16 train/valid/historical-test files to
  MACE, gates startup on archive preflight, and contains no blind-file argument
  or checkpoint-selection logic. Checkpoint review/selection remains with A.

## Findings for A

### Definite metadata defect: provenance counts sum to 29, not 35

`data/dataset_manifest.json` reports:

```text
split_sizes.train                         = 35
sum(provenance_status_counts_train)       = 29
sum(evidence_tier_counts_train)           = 29
new verified V16 archive additions        = 6
```

Both count maps are inherited from V15 and omit the six `verified_V16_PW_PBE`
rows. `preflight.py` validates split counts, finite labels, hashes, rank, and
holdout isolation, but does not check either provenance count map, so this
inconsistency does not block training preflight.

Exact reproduction:

```bash
python3 - <<'PY'
import json
p='research/mace-v12-transfer/periodic_interface_v4/mace_periodic_v16_interface_energy/data/dataset_manifest.json'
m=json.load(open(p))
print(m['split_sizes']['train'], sum(m['provenance_status_counts_train'].values()),
      sum(m['evidence_tier_counts_train'].values()))
PY
```

Recommended fix before treating these summaries as current: update the V16
builder to merge the six verified additions into both maps (adding a six-row
verified-V16 category), then make preflight assert both totals equal 35. If the
maps are intentionally parent-only, rename them accordingly and add separate
current-V16 totals. Keep the 34-row inherited `per_frame_provenance` ledger
scoped to V15; do not fabricate original-run evidence for the six newer rows.

### Nonblocking entry-point text mismatch

`run_v16_training.sh` is the V16 entry, but its usage line still names
`run_v15_training.sh` and its output-directory error message says “V15 run
directory”. `evaluate_v16.py` likewise labels its V16 output child as “V15 run
child”. The guards point at the V16 directory, so this is confusing text rather
than a path-selection failure. Update those strings to V16 before A hands the
entry point to another operator.

### Reproducibility gap: runtime policy scripts are not pinned

`input_pins.json` pins the data, build script, README, foundation model, and
frozen holdout inputs, but does not pin `run_v16_training.sh`, `preflight.py`,
or `evaluate_v16.py`. Their current contents were reviewed here, but a later
edit would not be detected by this pin set. Add hashes for those entry points
to the next reviewed V16 pin update.

## Provenance limits

The 29 inherited training rows each have one stored method declaration:
15 `source_method` fields and 14 `dft_method` fields, with no row carrying both
and none carrying neither. The parent provenance ledger still classifies 14 as
partial recovery, 12 as inherited declarations without original-run
reverification, and 3 as verified residual archives. The manifest correctly
keeps overall inherited physical-method consistency `UNKNOWN`. These gaps do
not invalidate finite serialized labels, but they prevent a claim that all
legacy labels have independently verified identical run settings.

## Artifacts

- `audit_v16_dataset.py`: reproducible read-only audit. The blind path is
  explicitly restricted to `pbe_interface_v16_blind_holdouts/inputs/`.
- `audit.json`: machine-readable checks, source hashes, method records, and
  identified metadata gaps. It contains no blind reference values.
- `SHA256SUMS.txt`: hashes for this report package.
