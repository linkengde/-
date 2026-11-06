#!/usr/bin/env bash
set -euo pipefail
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project="$(cd -- "$here/.." && pwd)"
binary="$project/tests/.p1p_kernel_test"
flat_binary="$project/tests/.p1p_flat_interface_test"
ag_binary="$project/tests/.p1p_ag_source_test"
trap 'rm -f -- "$binary" "$flat_binary" "$ag_binary"' EXIT
python3 "$project/tools/build_ag_evaporation_properties.py"
cc -std=c99 -O2 -Wall -Wextra -Werror "$here/test_p1p_evaporation_kernel.c" -lm -o "$binary"
cc -std=c99 -O2 -Wall -Wextra -Werror "$here/p1p_flat_interface_coupling_test.c" -lm -o "$flat_binary"
cc -std=c99 -O2 -Wall -Wextra -Werror "$here/test_ag_source_properties.c" -lm -o "$ag_binary"
cd "$project"
"$binary" "$here/synthetic_property_table.csv" "$project/tests/p1p_kernel_validation.csv"
"$flat_binary" "$here/synthetic_property_table.csv" "$project/tests/p1p_flat_interface_validation.csv"
"$ag_binary" "$project/data/Ag_evaporation_solver_table.csv" "$project/tests/p1p_ag_source_property_validation.csv"
