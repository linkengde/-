"""Non-scientific MPI test for the shared output-directory guard."""
import json
import os
from pathlib import Path

from gpaw.mpi import world

target = Path(os.environ["GPAW_MPI_GUARD_TEST_TARGET"])
target.mkdir(parents=True, exist_ok=False)
print(json.dumps({"rank": world.rank, "size": world.size, "target_exists": target.is_dir()}))
