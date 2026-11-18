#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="${HERE%/dft_kpoint_recheck_20261119}"
ACE_DIR="$HERE/ace_kconv_refits"
MTP_DIR="$HERE/mtp_kconv_refits"
OUT="$HERE/lammps_candidate_checks"
PACE_YAML2YACE="/workspace/.venvs/ace/bin/pace_yaml2yace"
MPIEXEC="/workspace/.venvs/qsgs/bin/mpiexec"
LMP_ACE="/workspace/potential-comparison/software/lammps-mpich-build/lmp"
LMP_MTP="/workspace/potential-comparison/software/lammps-mlip-build-src/src/lmp_mpi"
SMOKE_ACE="$BASE/ace_lammps_smoketest"
SMOKE_MTP="$BASE/mtp_lammps_smoketest"
mkdir -p "$OUT"

for weight in 1 0.05; do
  suffix="${weight/./}"
  model="$ACE_DIR/AgCSiTi_B1616_E10000_F1_a0.001_Ti_gap_weight${weight}.ace.yaml"
  dir="$OUT/ACE_weight_${suffix}"
  mkdir -p "$dir"
  cp "$SMOKE_ACE/data.lammps" "$dir/data.lammps"
  "$PACE_YAML2YACE" "$model" --output "$dir/model.yace" > "$dir/conversion.log" 2>&1
  cat > "$dir/in.ace" <<'EOF'
units metal
atom_style atomic
boundary p p p
read_data data.lammps
mass 1 107.8682
mass 2 12.011
mass 3 28.085
mass 4 47.867
pair_style pace
pair_coeff * * model.yace Ag C Si Ti
thermo 1
thermo_style custom step pe
run 0
write_dump all custom forces.dump id type x y z fx fy fz
EOF
  (cd "$dir" && "$MPIEXEC" -n 2 "$LMP_ACE" -in in.ace -log lammps.log > lammps.stdout 2>&1)
done

for tag in F005 F010 F020; do
  model="$MTP_DIR/MTP10_TiGap1x6x3_${tag}.mtp"
  [[ -f "$model" ]] || continue
  dir="$OUT/MTP_${tag}"
  mkdir -p "$dir"
  cp "$SMOKE_MTP/data.lammps" "$dir/data.lammps"
  cp "$model" "$dir/model.mtp"
  printf 'mtp-filename model.mtp\nselect FALSE\n' > "$dir/mlip.ini"
  cat > "$dir/in.mlip" <<'EOF'
units metal
atom_style atomic
boundary p p p
read_data data.lammps
mass 1 107.8682
mass 2 12.011
mass 3 28.085
mass 4 47.867
pair_style mlip mlip.ini
pair_coeff * *
thermo 1
thermo_style custom step pe
run 0
write_dump all custom forces.dump id type x y z fx fy fz
EOF
  (cd "$dir" && "$MPIEXEC" -n 2 "$LMP_MTP" -in in.mlip -log lammps.log > lammps.stdout 2>&1)
done
