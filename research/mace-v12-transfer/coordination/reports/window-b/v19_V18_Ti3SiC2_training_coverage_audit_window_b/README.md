# V18 training coverage for the Ti₃SiC₂ backbone

## Scope and lineage

This read-only audit covers the selected V18 clean-core `train.extxyz`. Before parsing it, the script verified its SHA-256 against both the selected-model record and dataset manifest: `20279d8aae53baf186546eb4a4dc396bb71e92de9d1c874d1c6fa0022313c045`. The selected epoch is 79. Frame provenance is taken from each training frame's own metadata. The only external structures read were the two public control geometries already screened by A: pure Ag and Ti₃SiC₂. Their source hashes and geometry identities match A's published screen.

No development, test, holdout, or sealed structure or label was opened. No other training dataset frames were opened. No MACE inference, training, DFT/MPI, data integration, label enabling, or edits to A-owned calculations occurred. The force inventory uses the `REF_forces` arrays already inside the selected training file.

## What is in the training set

There are 30 frames and 1,467 atom samples. Every frame contains Ag, Ti, Si, and C, and every frame is fully periodic. The formulas and configuration families are:

| Formula | Frames | Main family |
|---|---:|---|
| C13Ag17Si8Ti13 | 12 | Ag–Si |
| C13Ag2Si10Ti23 | 10 | Ag–C |
| C11Ag28SiTi9 | 7 | Ag–Ti |
| C8Ag10Si2Ti12 | 1 | Ti-gap |

There is **no pure Ti₃SiC₂ frame**, no pure Ag frame, and no frame whose non-Ag Ti:Si:C ratios are within 10% of 3:1:2. The closest family is C13Ag2Si10Ti23; its Ti/Si ratio is 23.3% below and its C/Si ratio is 35% below the ideal Ti₃SiC₂ ratio. The four species occur in these atom-sample counts across frames: Ag 430, Ti 461, Si 205, C 371.

Twenty-nine frames use large, approximately cubic 28 Å boxes; one is a compact, noncubic Ti-gap cell with lengths about 28.37, 3.08, and 5.33 Å. Repeated parent/source geometry hashes and paired displacement series show that 30 frames do not represent 30 independent structural families. All frames carry `proposal_role=train` and `proposal_provenance=directly_archive_verified`. Fifteen have an explicit source label/type, eight point to the same V14 training dataset path, and six carry parent geometry hashes that form three repeated pairs. Five frames also retain metadata saying they were scored as V13 geometry checks before entering V14 training; those rows are part of the current V18 training set and are not independent validation here. The per-frame source mapping is in `training_provenance.csv`; source archives outside the selected training file were not scanned.

## Local Ti–C environments

The A-screened Ti₃SiC₂ control is C16Si8Ti24 in a fully periodic hexagonal cell (6.115 × 6.115 × 17.6235 Å, γ=120°). Its Ti–C first shell has distances 2.088–2.176 Å; the next shell begins above 3.6 Å. We use 2.50 Å as the first-shell cutoff, inside that gap. The control has 16 Ti with three C neighbors, 8 Ti with six, and all 16 C with six Ti neighbors.

Training Ti–C distances overlap the first-shell range, but have a long tail: nearest-C distance per Ti has median 2.082 Å, p90 3.913 Å, and maximum 5.356 Å; nearest-Ti distance per C has median 2.065 Å, p90 4.495 Å, and maximum 6.418 Å. At the 2.50 Å cutoff across all training frames, Ti coordination counts are 62 atoms with zero C, 120 with one, 123 with two, 105 with three, 47 with four, and 4 with six. C coordination counts are 74 with zero Ti, 44 with one, 66 with two, 91 with three, 44 with four, 44 with five, and 8 with six.

The four fully six-coordinated Ti atoms and eight six-coordinated C atoms all come from the single compact Ti-gap frame. That frame is Ag-containing and chemically off-stoichiometric; it is not a pure Ti₃SiC₂ bulk example. None of the 29 Ag–C, Ag–Si, or Ag–Ti interface frames reaches six first-shell Ti–C neighbors on either species. Their local Ti/C coordination is dispersed and generally lower than the bulk control. The Ti-gap frame offers some similar local connectivity, while its composition and cell remain different from the screened pure phase.

The pure Ag control is Ag32 in a fully periodic 8.18 Å cubic cell; none of the training frames has that pure-phase geometry. A's existing screen reports force-vector RMSE 0.0000173 eV/Å for Ag and 2.88943 eV/Å for Ti₃SiC₂. The latter residual is almost entirely z-directed, with Ti and C species RMSEs 3.47676 and 2.62959 eV/Å. These are reported control metrics; this coverage audit did not rerun inference.

## Force components in `train.extxyz`

The existing training `REF_forces` contain substantial z components, so the file does not lack z-directed force magnitudes. Pooled over all frames:

| Species | Atom samples | x RMS (eV/Å) | y RMS (eV/Å) | z RMS (eV/Å) | max |Fz| (eV/Å) | Samples with |Fz| > 0.5 eV/Å |
|---|---:|---:|---:|---:|---:|---:|
| Ti | 461 | 1.708 | 1.640 | 1.441 | 5.647 | 304 |
| C | 371 | 1.769 | 1.417 | 1.501 | 4.014 | 271 |

The z-force values include both signs. They come from the mixed interface/contact configurations and do not establish matched z-displacement response for a chemically clean Ti₃SiC₂ backbone. Their magnitudes alone cannot explain or rule out the V18 error.

## Coverage diagnosis and proposed next labels

The data support a coverage-gap hypothesis: the V18 train set has many Ag-containing interface environments and Ti–C bonds at bulk-like distances, but lacks a clean, stoichiometric Ti₃SiC₂ training frame and lacks the bulk control's repeated Ti/C coordination pattern across the Ag–C/Si/Ti families. A single off-stoichiometric Ti-gap frame supplies some local six-coordination. Training labels also contain large z components, though not as a matched pristine-backbone response set. These counts describe coverage; they do not prove the microscopic cause of the model residual.

For a minimal follow-up acquisition proposal, without generating or registering structures:

1. Reuse the existing public Ti₃SiC₂ baseline geometry as the reference. Add four fixed-cell single points at the same PBE/PW 500 eV, σ=0.10 eV and 10×10×2 mesh: ±0.02 Å along z for one symmetry-distinct Ti and one C. This gives paired backbone response labels for the directions associated with the current force residual.
2. For Ag–Ti, Ag–Si and Ag–C, select one parent geometry per chemistry from a source structure family absent from the V18 training provenance. For each parent, evaluate the reference contact and ±0.10 Å normal separation, for nine single points. This is a first diagnostic panel; one parent per chemistry is not a statistical validation set, so further independent parents would be needed before a general validation claim.

No coordinates or DFT entries were generated, registered, or enabled. If resources require staging, prioritize the four Ti₃SiC₂ backbone points, then add the three interface families with independently sourced parents.

## Reproducibility

Run the read-only inventory from the repository root:

```bash
/workspace/.venvs/mace-cpu/bin/python research/mace-v12-transfer/coordination/reports/window-b/v19_V18_Ti3SiC2_training_coverage_audit_window_b/audit_train_coverage.py
```

The script checks hashes before it parses the training frames. Machine-readable summaries are in `inventory.json`; per-frame geometry and source fields, local coordination, and force distributions are in the CSV files. `SHA256SUMS.txt` covers every deliverable.
