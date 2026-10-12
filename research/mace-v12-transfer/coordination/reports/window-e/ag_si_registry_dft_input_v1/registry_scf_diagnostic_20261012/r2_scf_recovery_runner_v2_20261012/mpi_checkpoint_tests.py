#!/usr/bin/env python3
"""Four-rank synthetic checkpoint write-failure synchronization test."""
from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys
import tempfile

from ase import Atoms
from gpaw.mpi import broadcast, world

import runner


class RankZeroWriteFailure:
    def write(self, filename, mode=""):
        if mode != "all":
            raise AssertionError("wavefunctions were not requested")
        if world.rank == 0:
            raise OSError("synthetic rank-zero disk error")
        # GPAW's writer has a final collective barrier on successful peers.
        world.barrier()


def make_root() -> Path:
    name = tempfile.mkdtemp(prefix="e-r2-checkpoint-mpi-") if world.rank == 0 else None
    name = broadcast(name, 0, comm=world)
    return Path(name)


def main() -> int:
    if world.size != 4:
        raise SystemExit("test must be launched with exactly four MPI ranks")
    root = make_root()
    config, config_sha = runner.load_config()
    manager = runner.CheckpointManager(
        RankZeroWriteFailure(),
        Atoms("H", positions=[[0, 0, 0]], cell=[4, 4, 4], pbc=True),
        root, world, config, config_sha, runner.EXPECTED_LOCK_SHA,
        gpw_reader=lambda _path: {"wavefunctions": True},
    )
    caught = False
    try:
        manager.save(40)
    except runner.CheckpointFailure:
        caught = True
    caught_ranks = world.sum_scalar(int(caught))
    valid_artifact = (root / "checkpoint" / "checkpoint_iter_000040.gpw").exists()
    if world.rank == 0:
        result = {
            "test": "MPI4 rank-zero synthetic write error",
            "ranks": world.size,
            "all_ranks_stopped": caught_ranks == world.size,
            "valid_checkpoint_fabricated": valid_artifact,
            "status": "PASS" if caught_ranks == world.size and not valid_artifact else "FAIL",
        }
        (root / "write_failure.json").write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))
    world.barrier()
    if world.rank == 0:
        shutil.rmtree(root)
    world.barrier()
    return 0 if caught_ranks == world.size and not valid_artifact else 1


if __name__ == "__main__":
    raise SystemExit(main())
