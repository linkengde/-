#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
"/workspace/venvs/dft-labels/bin/python" "$root/sanitize_local_fit_data.py"
(cd "$root/POD" && "/workspace/software/lammps-pod-snap-clean/build/lmp" -in fitpod.in -log fitpod_v12g_forceorder_fixed.log)
test -s "$root/POD/pod_v12g_coefficients.pod"
test -s "$root/POD/pod_v12g_test_errors.pod"
(cd "$root/SNAP" && "/workspace/venvs/fitsnap/bin/python" /workspace/potential-fit/round7_pod_snap/run_fitsnap_py312.py --overwrite fitsnap.in > fitsnap_v12g.log 2>&1)
test -s "$root/SNAP/AgTiSiC_SNAP_v12g.snapcoeff"
test -s "$root/SNAP/AgTiSiC_SNAP_v12g.snapparam"
test -s "$root/SNAP/SNAP_v12g_metrics.md"
