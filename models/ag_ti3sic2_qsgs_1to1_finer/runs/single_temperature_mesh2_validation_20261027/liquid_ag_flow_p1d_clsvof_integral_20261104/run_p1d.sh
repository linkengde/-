#!/usr/bin/env bash
set -euo pipefail
: "${BASILISK_SRC:?Set BASILISK_SRC to the pinned Basilisk src directory}"
export BASILISK="$BASILISK_SRC"
export OPENGLIBS=-lfb_tiny
"$BASILISK_SRC/qcc" -O2 -Wall -autolink results/marangoni_Ca0066.c -o results/marangoni_Ca0066 -lm
timeout 900 results/marangoni_Ca0066 > results/Ca0066_stdout.log 2> results/Ca0066_stderr.log
"$BASILISK_SRC/qcc" -O2 -Wall -autolink results/marangoni_temperature_sigma_R16.c -o results/marangoni_temperature_sigma_R16 -lm
timeout 180 results/marangoni_temperature_sigma_R16 > results/temperature_sigma_R16.stdout.log 2> results/temperature_sigma_R16.stderr.log
