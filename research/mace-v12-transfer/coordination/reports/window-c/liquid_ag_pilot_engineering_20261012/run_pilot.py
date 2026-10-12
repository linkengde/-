#!/usr/bin/env python3
"""D's rho100 engineering-preflight CLI.

The fixed C lock is an execution HOLD. This command can run temporary
engineering fixtures but cannot turn the current lock into MD authorization.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from pilot_engine import (
    LOCK_SHA256,
    PreflightError,
    engineering_preflight,
    run_authorized,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("engineering-preflight", "run"), required=True)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--lock-sha256", default=LOCK_SHA256)
    parser.add_argument(
        "--potential-path",
        type=Path,
        default=Path("/workspace/WMH2006_Ag-local-validation/Ag.eam.alloy"),
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--authorization", type=Path)
    parser.add_argument("--authorization-sha256")
    args = parser.parse_args()
    try:
        if args.mode == "engineering-preflight":
            result = engineering_preflight(args.lock, args.potential_path)
            output = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
            if args.output_json:
                args.output_json.parent.mkdir(parents=True, exist_ok=True)
                args.output_json.write_text(output, encoding="utf-8")
            sys.stdout.write(output)
            return 0
        if args.output_dir is None:
            raise PreflightError("--output-dir is required for run mode")
        if args.authorization is None or not args.authorization_sha256:
            raise PreflightError(
                "run mode requires C's separately published execution authorization and SHA-256"
            )
        result = run_authorized(
            args.lock, args.lock_sha256, args.potential_path, args.output_dir,
            args.authorization, args.authorization_sha256,
        )
        sys.stdout.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
        return 0
    except PreflightError as exc:
        output = {
            "PREFLIGHT_STATUS": "FAIL",
            "execution_authorized": False,
            "ERROR": str(exc),
        }
        sys.stderr.write(json.dumps(output, indent=2, sort_keys=True) + "\n")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
