"""Rank-0-only directory creation guard for the frozen GPAW input script."""
import os
import time
from pathlib import Path

_original_mkdir = Path.mkdir
_enabled = os.environ.get("GPAW_MPI_OUTPUT_CREATE_GUARD") == "1"
_rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", os.environ.get("PMI_RANK", os.environ.get("PMIX_RANK", "0"))))


def _mkdir(self, mode=0o777, parents=False, exist_ok=False):
    if _enabled and self.name == "coarse" and self.parent.name == "results" and not exist_ok:
        if _rank == 0:
            return _original_mkdir(self, mode=mode, parents=parents, exist_ok=exist_ok)
        deadline = time.monotonic() + 60.0
        while not self.exists():
            if time.monotonic() >= deadline:
                raise TimeoutError(f"rank {_rank} timed out waiting for rank 0 to create {self}")
            time.sleep(0.02)
        return None
    return _original_mkdir(self, mode=mode, parents=parents, exist_ok=exist_ok)


if _enabled:
    Path.mkdir = _mkdir
