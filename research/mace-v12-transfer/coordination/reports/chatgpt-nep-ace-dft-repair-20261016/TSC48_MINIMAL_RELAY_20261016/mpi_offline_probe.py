"""No calculator/DFT: exercise four-rank adapter collectives with synthetic bytes only."""
import sys,time
from pathlib import Path
from gpaw.mpi import world
import relay_runner as relay
folder=Path(sys.argv[1])
class FakeCalculator:
    def write(self,path,mode):
        assert mode=='all'
        if world.rank==0:Path(path).write_bytes(b'SYNTHETIC_NOT_GPAW'*1024)
        world.barrier()
for n in (1,2):
    relay.mpi_save(FakeCalculator(),world,folder,{'scope':'MPI_SYNTHETIC_ONLY'},n,n,
                   None,n,time.monotonic(),0,n,'synthetic.log')
try:
    relay.mpi_save(FakeCalculator(),world,folder,{'scope':'MPI_SYNTHETIC_ONLY'},3,3,
                   None,3,time.monotonic(),0,3,'synthetic.log')
except relay.base.Blocked:
    print('rank=%d TWO_FILES_AND_THIRD_WRITE_BLOCKED=PASS'%world.rank)
else:raise AssertionError('third write should be blocked')
