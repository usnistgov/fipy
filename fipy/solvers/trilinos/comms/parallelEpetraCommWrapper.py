from __future__ import unicode_literals
from mpi4py import MPI

from fipy.tools import numerix
from fipy.solvers.trilinos.comms.epetraCommWrapper import EpetraCommWrapper

__all__ = ["ParallelEpetraCommWrapper"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class ParallelEpetraCommWrapper(EpetraCommWrapper):
    """MPI Communicator wrapper

    Encapsulates capabilities needed for both Epetra and mpi4py.
    """
    
    def __init__(self):
        self.mpi4py_comm = MPI.COMM_WORLD
        super(ParallelEpetraCommWrapper, self).__init__()
        
    def __setstate__(self, dict):
        self.__init__()
        
    def all(self, a, axis=None):
        return self.mpi4py_comm.allreduce(a.all(axis=axis), op=MPI.LAND)

    def any(self, a, axis=None):
        return self.mpi4py_comm.allreduce(a.any(axis=axis), op=MPI.LOR)

    def allclose(self, a, b, rtol=1.e-5, atol=1.e-8):
        return self.mpi4py_comm.allreduce(numerix.allclose(a, b, rtol=rtol, atol=atol), op=MPI.LAND)

    def allequal(self, a, b):
        return self.mpi4py_comm.allreduce(numerix.allequal(a, b), op=MPI.LAND)

    def bcast(self, obj, root=0):
        return self.mpi4py_comm.bcast(obj=obj, root=root)

    def allgather(self, obj):
        """mpi4py `allgather`
        
        Communicates copies of each `sendobj` to every rank in the comm, creating
        a rank-dimensional list of `sendobj` objects.
        
        >>> m4count = self.mpi4py_comm.allgather(self.mpi4py_comm.Get_rank())
        >>> from builtins import range
        >>> for i in range(self.mpi4py_comm.Get_size()):
        ...     assert m4count[i] == i

        """
        return self.mpi4py_comm.allgather(sendobj=obj)

    def MaxAll(self, obj):
        """return max across all processes
        """
        return self.mpi4py_comm.allreduce(sendobj=obj, op=max)

    def MinAll(self, obj):
        """return min across all processes
        """
        return self.mpi4py_comm.allreduce(sendobj=obj, op=min)
