"""Read-only check of GPAW 26.7.0's remaining-SCF-iteration API."""
import argparse
import json

from gpaw import restart
from gpaw.mpi import world

parser = argparse.ArgumentParser()
parser.add_argument("--checkpoint", required=True)
parser.add_argument("--prior-iterations", type=int, required=True)
parser.add_argument("--locked-maxiter", type=int, required=True)
args = parser.parse_args()

atoms, calc = restart(args.checkpoint, txt=None)
before = int(calc.dft.scf_loop.maxiter)
remaining = args.locked_maxiter - args.prior_iterations
if before != args.locked_maxiter or remaining <= 0:
    raise RuntimeError(f"Unexpected checkpoint loop state: {before=}, {remaining=}")
calc.dft.scf_loop.maxiter = remaining
after = int(calc.dft.scf_loop.maxiter)
if after != remaining:
    raise RuntimeError(f"Failed to assign remaining iteration cap: {after=}")

if world.rank == 0:
    print(json.dumps({
        "probe": "PASS_NO_SCF",
        "atom_count": len(atoms),
        "checkpoint_iteration": args.prior_iterations,
        "checkpoint_stored_maxiter": before,
        "locked_total_maxiter": args.locked_maxiter,
        "remaining_iterations_assigned": after,
        "scf_started": False,
    }, sort_keys=True))
world.barrier()
