#!/usr/bin/env bash
set -euo pipefail
case_dir=$(pwd)
[[ -f run_config.json && -f latest_checkpoint.txt ]] || { echo 'Run from a case directory with run_config.json and latest_checkpoint.txt' >&2; exit 2; }
for p in /proc/[0-9]*; do
  [[ -e "$p/cwd" && -r "$p/cmdline" ]] || continue
  cwd=$(readlink "$p/cwd" 2>/dev/null || true)
  [[ "$cwd" == "$case_dir" ]] || continue
  cmd=$(tr '\0' ' ' < "$p/cmdline" 2>/dev/null || true)
  if [[ "$cmd" == *p1g_free_surface_coupled_resume_safe* ]]; then
    echo "A P1G solver is already running in $case_dir (PID ${p##*/}); refusing duplicate resume." >&2
    exit 3
  fi
done
checkpoint=$(<latest_checkpoint.txt)
if [[ ! -f "$checkpoint" && -f "$checkpoint.gz" ]]; then
  pigz -dc -- "$checkpoint.gz" > "$checkpoint.part"
  mv "$checkpoint.part" "$checkpoint"
fi
[[ -f "$checkpoint" && -f "$checkpoint.state" && -f "$checkpoint.complete" ]] || {
  echo "Checkpoint is incomplete: $checkpoint" >&2; exit 4;
}
mode=$(python3 -c 'import json; print(json.load(open("run_config.json"))["marangoni"])')
[[ "$mode" == on || "$mode" == off ]] || { echo "Invalid Marangoni mode: $mode" >&2; exit 5; }
if [[ -f exit_code.txt ]]; then
  mv exit_code.txt "prior_exit_code_$(date -u +%Y%m%dT%H%M%SZ).txt"
fi
export P1G_RESTART_DUMP="$case_dir/$checkpoint"
../../watch_p1g_outputs.sh >"compression_watch_resume_$(date -u +%Y%m%dT%H%M%SZ).log" 2>&1 &
compressor_pid=$!
date -u '+%Y-%m-%dT%H:%M:%SZ' >> resume_started_utc.txt
set +e
timeout --signal=TERM --kill-after=60s 259200s env P1G_MARANGONI="$mode" ../../p1g_free_surface_coupled_resume_safe ../../../inputs/basilisk negative 1e11 1000000 >>run.stdout.log 2>>run.stderr.log
status=$?
set -e
printf '%s\n' "$status" > exit_code.txt
wait "$compressor_pid" || true
exit "$status"
