


__docformat__ = 'restructuredtext'

from fipy.tools import numerix

__all__ = ["CommWrapper", "ParallelCommWrapper"]

class CommWrapper(object):
    """MPI Communicator wrapper

    Encapsulates capabilities needed for Epetra. Some capabilities are not parallel.

    """

    def __init__(self, Epetra=None):
        self.epetra_comm = Epetra.PyComm()

    def __repr__(self):
        return "%s()" % self.__class__.__name__

    @property
    def procID(self):
        return self.epetra_comm.MyPID()

    @property
    def Nproc(self):
        return self.epetra_comm.NumProc()

    def Barrier(self):
        self.epetra_comm.Barrier()

    def all(self, a, axis=None):
        return a.all(axis=axis)

    def any(self, a, axis=None):
        return a.any(axis=axis)

    def allclose(self, a, b, rtol=1.e-5, atol=1.e-8):
        return numerix.allclose(a, b, rtol=rtol, atol=atol)

    def allequal(self, a, b):
        return numerix.allequal(a, b)

    def bcast(self, obj, root=0):
        return obj

    def allgather(self, obj):
        return obj

    def sum(self, a, axis=None):
        summed = numerix.array(a).sum(axis=axis)
        shape = summed.shape
        if shape == ():
            summed = summed.reshape((1,))
        parallelSummed = self.epetra_comm.SumAll(summed)
        if shape == ():
            parallelSummed = parallelSummed.reshape(())
        return parallelSummed

    def __getstate__(self):
        return {'dummy': 0}

    def __setstate__(self, dict):
        from PyTrilinos import Epetra
        self.__init__(Epetra=Epetra)

    def Norm2(self, vec):
        return vec.Norm2()

    def MaxAll(self, vec):
        return self.epetra_comm.MaxAll(numerix.array(vec))

    def MinAll(self, vec):
        return self.epetra_comm.MinAll(numerix.array(vec))

class ParallelCommWrapper(CommWrapper):
    """MPI Communicator wrapper for parallel processes"""
    pass
