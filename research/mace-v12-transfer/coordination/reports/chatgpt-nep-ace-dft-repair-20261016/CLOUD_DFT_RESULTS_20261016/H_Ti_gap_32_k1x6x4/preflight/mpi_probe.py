from gpaw.mpi import world
assert world.size == 4
if world.rank == 0:
    print('GPAW_MPI_4_RANK_IMPORT_PASS', flush=True)
