# Result audit: TSC96_DEV_PRISM_10m10_COD9009647_8L_RESP_TIZ020

Status: **SUCCESS**. The fixed-geometry GPAW single point converged in 52 SCF iterations. `RUN_STATUS.json`, `RESULT.json`, the GPAW log, and ASE-read labels agree on convergence and saved quantities.

## Binding and geometry checks

- Input SHA-256: `710b0abdd20f0d59103ebfe061ad65c834bae7a4b2f5a300453ba1527e63e030` (matches the frozen input file).
- Geometry SHA-256: `33aa1032f0ec404b5ca46976f82da5c4e1dc29cff8cebb76cf3f25d8a7a92e49`.
- Scientific configuration SHA-256: `c1b19b90a091d735bbbc6ebb184fb31ec1f7d3985c75f952879bb43bfda85fba`.
- Runner SHA-256: `93d0b72fca26fadbdc3b6006dcc65341a0452b4c9672f6f2e6a468981a27c7a9`.
- Atoms: 96 (Ti48 Si16 C32); symbols/order, positions, cell, and PBC in `labelled.extxyz` are exactly identical to the input when read through ASE.
- Fixed coordinates and cell; no relaxation was performed.
- All 96 force vectors, energy, force-consistent free energy, and the 3×3 stress tensor are finite. Values parsed from `labelled.extxyz` match `RESULT.json` (force serialization agreement within 1e-7 eV/Å; energy/free energy and stress within the checks recorded by the validator).
- GPAW 26.7.0; ASE 3.29.0; PBE/PAW, PW 500 eV; Fermi–Dirac width 0.10 eV; standard Monkhorst–Pack 12×2×1 with no explicit Gamma flag or added shift; fixed geometry. Four MPI ranks, one OpenMP thread, from the launch preflight.

## Results

- Total energy: **-752.548144181147 eV**.
- Force-consistent free energy: **-753.272762694337 eV**.
- Maximum atomic force magnitude: **2.01217427 eV/Å**.
- Stress (eV/Å³): `[[-0.0027331425426275428, -1.5787404115500424e-06, 2.0952601502293885e-09], [-1.5787404115500424e-06, -0.004060401298185841, -4.275033873913782e-05], [2.0952601502293885e-09, -4.275033873913782e-05, -0.0008650592367078775]]`.
- Wall time from GPAW timing: **16,503.036 s** (about 4.58 h); SCF portion 15,857.930 s.
- Peak sampled cgroup memory from the resource profile: **16589258752 bytes** (includes file cache in the cgroup accounting).

## Data-role limitation

Candidate role remains `DEV_GEOMETRY_CANDIDATE_TITANIUM_DISPLACEMENT_RESPONSE`; no TRAIN/DEV/TEST assignment was changed. This displacement structure belongs to source family `surface_prismatic_10m10_cod9009647_8layer` and shares COD 9009647 with the existing reference and related TRAIN source. It can support a geometry-response candidate analysis, but it is not an independent crystallographic-source blind test.

## Temporary checkpoint

The successful result and labels were validated before checkpoint disposal. The temporary checkpoint had size 2422680382 bytes and SHA-256 `cf086421df3e017c58f960639c220db6063ab1f912602da23ac0c0657d082b61`. It was not resumed or needed for recovery; its disposition is recorded separately. The checkpoint file itself is intentionally excluded from the archival package.
