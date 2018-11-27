


from fipy.tools.comms.commWrapper import CommWrapper
from fipy.tools import numerix

__all__ = ["Mpi4pyCommWrapper"]

class Mpi4pyCommWrapper(CommWrapper):
    """MPI Communicator wrapper

    Encapsulates capabilities needed for both Epetra and mpi4py.

    """

    def __init__(self, Epetra, MPI):
        self.MPI = MPI
        self.mpi4py_comm = self.MPI.COMM_WORLD
        CommWrapper.__init__(self, Epetra)

    def __setstate__(self, dict):
        from PyTrilinos import Epetra
        from mpi4py import MPI
        self.__init__(Epetra=Epetra, MPI=MPI)

    def all(self, a, axis=None):
        return self.mpi4py_comm.allreduce(a.all(axis=axis), op=self.MPI.LAND)

    def any(self, a, axis=None):
        return self.mpi4py_comm.allreduce(a.any(axis=axis), op=self.MPI.LOR)

    def allclose(self, a, b, rtol=1.e-5, atol=1.e-8):
        return self.mpi4py_comm.allreduce(numerix.allclose(a, b, rtol=rtol, atol=atol), op=self.MPI.LAND)

    def allequal(self, a, b):
        return self.mpi4py_comm.allreduce(numerix.allequal(a, b), op=self.MPI.LAND)

    def bcast(self, obj, root=0):
        return self.mpi4py_comm.bcast(obj=obj, root=root)

    def allgather(self, obj):
        """mpi4py allgather
        
        Communicates copies of each sendobj to every rank in the comm, creating
        a rank-dimensional list of sendobj objects.
        
        >>> m4count = self.mpi4py_comm.allgather(self.mpi4py_comm.Get_rank())
        >>> for i in range(self.mpi4py_comm.Get_size()):
        ...     assert m4count[i] == i

        """
        return self.mpi4py_comm.allgather(sendobj=obj)
