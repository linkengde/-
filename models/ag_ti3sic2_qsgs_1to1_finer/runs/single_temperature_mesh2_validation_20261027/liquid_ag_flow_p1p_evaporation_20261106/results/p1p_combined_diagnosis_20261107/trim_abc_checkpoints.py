#!/usr/bin/env python3
"""Keep the two newest verified checkpoints in one isolated A/B/C run folder."""
import re
import sys
from pathlib import Path

folder = Path(sys.argv[1])
evaporation_enabled = len(sys.argv) > 2 and sys.argv[2] == "1"
pattern = re.compile(r"^p1h_restart_i(\d+)_t([0-9.]+)_ps\.dump$")
valid = []
for dump in folder.glob("p1h_restart_i*_t*_ps.dump"):
    m = pattern.match(dump.name)
    if not m:
        continue
    bundle = [dump, Path(str(dump) + ".complete"), Path(str(dump) + ".state")]
    p1p = Path(str(dump) + ".p1p")
    if evaporation_enabled:
        bundle.append(p1p)
    if all(p.is_file() and p.stat().st_size for p in bundle):
        valid.append((float(m.group(2)), int(m.group(1)), bundle))
valid.sort()
for _, _, bundle in valid[:-2]:
    for path in bundle:
        path.unlink(missing_ok=True)
