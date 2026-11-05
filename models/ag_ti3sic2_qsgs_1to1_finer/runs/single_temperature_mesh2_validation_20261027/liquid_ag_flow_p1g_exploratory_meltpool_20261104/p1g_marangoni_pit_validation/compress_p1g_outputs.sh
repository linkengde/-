#!/usr/bin/env bash
set -euo pipefail
latest=""
if [[ -f latest_checkpoint.txt ]]; then latest=$(<latest_checkpoint.txt); fi
for f in p1g_*_t*_ps.bin; do
  [[ -f "$f" ]] || continue
  pigz -1 -n -- "$f"
done
for f in p1g_restart_*.dump; do
  [[ -f "$f" ]] || continue
  [[ -f "$f.complete" ]] || continue
  [[ "$f" == "$latest" ]] && continue
  pigz -1 -n -- "$f"
done
