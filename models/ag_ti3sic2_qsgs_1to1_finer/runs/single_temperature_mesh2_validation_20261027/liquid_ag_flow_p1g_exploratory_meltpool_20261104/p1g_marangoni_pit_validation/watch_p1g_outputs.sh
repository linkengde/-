#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
while [[ ! -f exit_code.txt ]]; do
  "$script_dir/compress_p1g_outputs.sh"
  sleep 10
done
"$script_dir/compress_p1g_outputs.sh"
