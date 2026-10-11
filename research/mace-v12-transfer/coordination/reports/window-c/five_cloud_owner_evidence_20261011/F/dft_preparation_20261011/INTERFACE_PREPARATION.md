# Ti₃SiC₂/Ag interface DFT preparation

This document prepares an interface input record only. No interface geometry is supplied or generated here.

## Assignment boundary

F's six task IDs in `HANDOFF_F.md` are liquid Ag candidates, pure Ti₃SiC₂ response, and solid Ag response. None is a mixed Ti₃SiC₂/Ag interface task. The existing C plan assigns future Ag-framework interface work to a separate B2 proposal/owner flow. Before preparing a runnable mixed-interface input, C must assign an interface task ID and approve its parent geometry and role. Do not attach an interface calculation to `C-F-TSC-PF01-CFG03` or `C-F-AG-PF01-CFG03`.

## Geometry evidence required from the owner

The interface manifest must identify and hash the source structures and all transformations, then freeze:

- Ti₃SiC₂ facet/orientation, termination and source cell; Ag orientation/source cell;
- in-plane supercell match, integer transformation matrices, residual strain on each side, and the rule used to distribute strain;
- interface registry, lateral shift, Ag-framework separation definition and its source, slab thickness, atom count, ordered atom IDs, and element order;
- cell vectors, PBC, vacuum thickness, surface normal, fixed-layer mask (if any), and any dipole correction axis;
- parent-family ID, source CIF/structure hashes, provenance/license review, geometry hash, and role/split assignment.

Do not select termination, registry, gap, strain, vacuum, or k points by inspecting new DFT labels. The C plan's suggested Ti-terminated basal, Si-terminated basal and mixed prismatic motifs are design candidates only; they are not frozen F inputs. Its `D0` separation must come from a provenance-reviewed source before any gap offsets are defined.

## Method preparation

For a later GPAW interface calculation, preserve the project baseline: PBE, PAW data 1.2.1, PW 500 eV, Fermi-Dirac sigma 0.10 eV, and the reviewed SCF thresholds. Set k points from the frozen in-plane cell and a reviewed convergence comparison; use one point along the non-periodic axis only after that geometry and PBC are approved. The current project plan calls for at least 20–26 Å slab vacuum and same-geometry vacuum/k-point checks. Record the actual dipole treatment and its direction, not just a Boolean flag.

An interface has xy periodicity and vacuum along z only if the approved cell says so. Do not infer PBC, surface normal or dipole axis from filenames. Do not run cell relaxation or ionic relaxation under a static single-point task.

## Outputs and interpretation

Save the complete per-atom forces and the raw symmetric 3×3 ASE stress tensor in eV/Å³ with its Voigt representation `[xx, yy, zz, yz, xz, xy]`. For a slab with vacuum, label that as raw cell stress; do not interpret the vacuum-diluted three-dimensional stress as a bulk stress. Surface stress requires a separately defined area normalization and method.

An interface total energy alone is not an interface formation energy. Any later interface-energy comparison needs separately reviewed, method-matched reference calculations, consistent composition/atom counts, and an explicit normalization (for example per interface area). No interface energy is computed in this preparation.

The placeholder fields are in `INTERFACE_INPUT_MANIFEST.template.json`.
