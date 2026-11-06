#!/usr/bin/env bash
set -euo pipefail
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project="$(cd -- "$here/.." && pwd)"
binary="$project/tests/.p1p_kernel_test"
flat_binary="$project/tests/.p1p_flat_interface_test"
trap 'rm -f -- "$binary" "$flat_binary"' EXIT
cc -std=c99 -O2 -Wall -Wextra -Werror "$here/test_p1p_evaporation_kernel.c" -lm -o "$binary"
cc -std=c99 -O2 -Wall -Wextra -Werror "$here/p1p_flat_interface_coupling_test.c" -lm -o "$flat_binary"
cd "$project"
"$binary" "$here/synthetic_property_table.csv" "$project/tests/p1p_kernel_validation.csv"
"$flat_binary" "$here/synthetic_property_table.csv" "$project/tests/p1p_flat_interface_validation.csv"
