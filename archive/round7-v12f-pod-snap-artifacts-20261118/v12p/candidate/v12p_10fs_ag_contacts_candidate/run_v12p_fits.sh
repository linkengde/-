#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
(cd "$root/POD" && "/workspace/software/lammps-pod-snap-clean/build/lmp" -in fitpod.in -log fitpod_v12p.log)
test -s "$root/POD/pod_v12p_coefficients.pod"
test -s "$root/POD/pod_v12p_test_errors.pod"
(cd "$root/SNAP" && "/workspace/venvs/fitsnap/bin/python" /workspace/potential-fit/round7_pod_snap/run_fitsnap_py312.py --overwrite fitsnap.in > fitsnap_v12p.log 2>&1)
test -s "$root/SNAP/AgTiSiC_SNAP_v12p.snapcoeff"
test -s "$root/SNAP/AgTiSiC_SNAP_v12p.snapparam"
test -s "$root/SNAP/SNAP_v12p_metrics.md"
