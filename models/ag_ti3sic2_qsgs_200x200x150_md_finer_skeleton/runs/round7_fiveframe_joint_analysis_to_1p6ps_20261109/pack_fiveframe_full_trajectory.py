#!/usr/bin/env python3
"""Pack five complete Round7 LAMMPS snapshots into one OVITO-playable dump.gz."""
from __future__ import annotations
import gzip
from pathlib import Path

MODEL = Path(__file__).resolve().parents[2]
OLD_RUN = MODEL / "runs" / "ttm_mtp_round7_continue_to_1p5ps_20261109"
NEW_RUN = MODEL / "runs" / "ttm_mtp_round7_continue_1p4_to_1p6ps_20261109"
OUT = Path(__file__).resolve().parent / "round7_5frames_full_system_0p6_to_1p6ps.dump.gz"
TIMES = (0.60, 0.85, 1.10, 1.35, 1.60)
EXPECTED_ATOMS = 352_637


def frame_path(time_ps: float) -> Path:
    tag = f"{time_ps:.2f}".rstrip("0").rstrip(".").replace(".", "p")
    runs = (NEW_RUN, OLD_RUN)
    for run in runs:
        path = run / "checkpoints" / f"atoms_{tag}ps.dump"
        if path.is_file():
            return path
    raise FileNotFoundError(f"No whole-system frame at {time_ps:.2f} ps")


def validate(path: Path, expected_step: int) -> None:
    with path.open("rt", encoding="ascii") as stream:
        header = [next(stream) for _ in range(9)]
        timestep = int(header[1].strip())
        count = int(header[3].strip())
        if timestep != expected_step or count != EXPECTED_ATOMS:
            raise ValueError(f"Unexpected snapshot {path}: step={timestep}, atoms={count}")
        lines = sum(1 for _ in stream)
        if lines != count:
            raise ValueError(f"Incomplete snapshot {path}: got {lines}/{count} atom rows")


def main() -> None:
    sources = [(time_ps, frame_path(time_ps)) for time_ps in TIMES]
    for time_ps, path in sources:
        validate(path, round(time_ps * 10_000))
        print(f"validated {time_ps:.2f} ps: {path}", flush=True)
    with gzip.open(OUT, "wt", encoding="ascii", compresslevel=3) as target:
        for _, path in sources:
            with path.open("rt", encoding="ascii") as source:
                while chunk := source.read(4 * 1024 * 1024):
                    target.write(chunk)
    print(f"wrote {OUT} ({OUT.stat().st_size} bytes)", flush=True)


if __name__ == "__main__":
    main()
