#!/usr/bin/env python3
"""Read-only R2 recovery runtime, source, CPU, memory and storage preflight."""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import runner  # noqa: E402


def main() -> int:
    payload = runner.run_preflight(Path("/workspace"))
    payload["checked_at_utc"] = datetime.now(timezone.utc).isoformat()
    payload["environment"] = "Codex Cloud task-isolated temporary workspace"
    payload["current_git_branch"] = "cloud-e-b1-dft"
    payload["current_git_commit"] = "19d9b672cc75c6297646d01f491427cf581b2569"
    output = HERE / "RESOURCE_PREFLIGHT.json"
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0 if payload["status"] == "PREFLIGHT_PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
