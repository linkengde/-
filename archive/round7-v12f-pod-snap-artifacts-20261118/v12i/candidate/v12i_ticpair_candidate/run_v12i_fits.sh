#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
(cd "$root/POD" && "/workspace/software/lammps-pod-snap-clean/build/lmp" -in fitpod.in -log fitpod_v12i.log)
test -s "$root/POD/pod_v12i_coefficients.pod"
test -s "$root/POD/pod_v12i_test_errors.pod"
(cd "$root/SNAP" && "/workspace/venvs/fitsnap/bin/python" /workspace/potential-fit/round7_pod_snap/run_fitsnap_py312.py --overwrite fitsnap.in > fitsnap_v12i.log 2>&1)
test -s "$root/SNAP/AgTiSiC_SNAP_v12i.snapcoeff"
test -s "$root/SNAP/AgTiSiC_SNAP_v12i.snapparam"
test -s "$root/SNAP/SNAP_v12i_metrics.md"
