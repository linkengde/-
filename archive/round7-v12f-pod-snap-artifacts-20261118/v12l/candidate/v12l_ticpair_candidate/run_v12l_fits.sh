#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
(cd "$root/POD" && "/workspace/software/lammps-pod-snap-clean/build/lmp" -in fitpod.in -log fitpod_v12l.log)
test -s "$root/POD/pod_v12l_coefficients.pod"
test -s "$root/POD/pod_v12l_test_errors.pod"
(cd "$root/SNAP" && "/workspace/venvs/fitsnap/bin/python" /workspace/potential-fit/round7_pod_snap/run_fitsnap_py312.py --overwrite fitsnap.in > fitsnap_v12l.log 2>&1)
test -s "$root/SNAP/AgTiSiC_SNAP_v12l.snapcoeff"
test -s "$root/SNAP/AgTiSiC_SNAP_v12l.snapparam"
test -s "$root/SNAP/SNAP_v12l_metrics.md"
