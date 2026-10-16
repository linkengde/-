#!/usr/bin/env bash
# Run in a genuinely independent H2 cloud window; synthetic data only, no DFT.
set -euo pipefail
relay_test_dir=$(mktemp -d /tmp/tsc48-h2-readback.XXXXXX)
relay_base='https://raw.githubusercontent.com/linkengde/-/bdf0f292edc1d75bbfe41655f70373838ca9a4ab/research/mace-v12-transfer/coordination/reports/chatgpt-nep-ace-dft-repair-20261016/TSC48_MINIMAL_RELAY_20261016'
for relay_file in chunk_transfer.py stream_join_github.py; do
  curl -fsSL "$relay_base/$relay_file" -o "$relay_test_dir/$relay_file"
done
curl -fsSL "$relay_base/synthetic_parts/manifest.json" -o "$relay_test_dir/manifest.json"
printf '%s  %s\n' b990e993e246f14c6762e01754f685c140d7a26e3333d1fd069806815bab49e6 "$relay_test_dir/manifest.json" | sha256sum -c -
python "$relay_test_dir/stream_join_github.py" "$relay_test_dir/manifest.json" "$relay_base/synthetic_parts/" "$relay_test_dir/restored.bin"
printf '%s  %s\n' ce41df6ec4ca386de691755fa545f4b2f9255b95b186795c8bf40910cdd79022 "$relay_test_dir/restored.bin" | sha256sum -c -
printf 'H2_OUTPUT_DIRECTORY=%s\n' "$relay_test_dir"
