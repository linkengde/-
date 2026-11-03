#!/usr/bin/env bash
set -euo pipefail

pilot_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
case_name="${1:-static_drop_case}"
solver="${2:-interFoam}"
case_dir="$pilot_dir/$case_name"
image="opencfd/openfoam-default:2406"

if [[ ! -d "$case_dir" ]]; then
    printf 'case directory not found: %s\n' "$case_dir" >&2
    exit 2
fi

if [[ "$solver" != "interFoam" && "$solver" != "interIsoFoam" ]]; then
    printf 'unsupported solver: %s\n' "$solver" >&2
    exit 2
fi

docker --host=unix:///var/run/docker.sock run --rm \
    --user "$(id -u):$(id -g)" \
    --mount "type=bind,src=$case_dir,dst=/case" \
    --workdir /case \
    --entrypoint bash "$image" -lc \
    'blockMesh > log.blockMesh 2>&1 && "$1" > "log.$1" 2>&1' bash "$solver"
