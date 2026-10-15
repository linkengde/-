# C review request: Ti₃SiC₂ P0 DFT response batch

The NEP/ACE data audit found that the existing 48-atom Ti₃SiC₂ parent has no labelled volume/shear stress response or Si displacement pair. E has prepared six geometry-only single-point candidates in `ti3sic2_bulk_p0/`; no calculation has been run.

## Candidate batch

- `TSC48_VOL_m02`, `TSC48_VOL_p02`: uniform volume factors 0.98 and 1.02.
- `TSC48_BASAL_XY_SHEAR_m005`, `TSC48_BASAL_XY_SHEAR_p005`: basal-plane simple shear `x′ = x + γy`, γ = −0.005/+0.005.
- `TSC48_SI_ID7_Z_m020`, `TSC48_SI_ID7_Z_p020`: local Si ID 7, z displacement −0.02/+0.02 Å.

All six derive from one source parent (48 atoms, Ti24Si8C16) and must remain together in any future data role. Existing `role_freeze.json` is unchanged: TRAIN/DEV/TEST remain empty. No role is assigned by this request.

## Proposed method (not locked)

The candidate lock proposes the common archived PBE/PAW setup: GPAW 26.7.0, ASE 3.29.0, 500 eV plane-wave cutoff, Fermi–Dirac 0.10 eV, SCF energy/density/eigenstate thresholds 1e-5, Pulay beta 0.05/nmaxold 8/weight 100, and maxiter 160. kpts 10×10×4 is the highest archived common bulk/displacement mesh. The archived input metadata does not establish k-point centering; MPI rank count, wall-time limit, checkpoint policy and work root are not frozen. These must be decided before any run.

Request output for fixed-geometry single points: native total energy, force-consistent free energy, per-atom forces, stress tensor for volume/shear points, SCF convergence and iteration count, plus input/result SHA-256. Keep each structure and cell fixed; no relaxation.

## Environment facts

The pre-existing E-private environment passes software/PAW identity checks: GPAW 26.7.0, ASE 3.29.0, gpaw-data 1.2.1, OpenMPI 5.0.7, MPI enabled; PAW hashes match the project manifest. E runtime configures `GPAW_SETUP_PATH` and `GPAW_MPI_BACKEND=cgpaw`. The current Codex workspace has 4 CPU quota cores, 32 GiB memory and about 22.43 GB free, but `/mnt/e-work` is not mounted and `DFT_WORK_ROOT` is unset. The environment is CPU-only.

## Decisions requested from C

1. Approve or amend these six exact structures and the basal-shear convention.
2. Freeze the role of the entire 48-atom parent family, or explicitly retain it as diagnostic-only; do not split its children across roles.
3. Freeze k-point centering, MPI/OMP, wall-time/checkpoint policy and output contract.
4. Select an approved DFT work root that fits this workspace's temporary-storage semantics.
5. If accepted, issue a separate execution authorization bound to the structure-manifest and parameter-proposal SHA values.

Until those decisions and authorization are present, E keeps the batch at `PREPARED_NOT_RUN`. Further surface/interface/liquid/contact structures remain separate priorities in the audit; these six points alone do not qualify NEP/ACE training.
