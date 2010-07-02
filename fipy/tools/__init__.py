try:
    import scipy
except:
    pass

try:
    from PyTrilinos import Epetra
    from mpi4py import MPI
    
    class CommWrapper(object):
        """MPI Communicator wrapper
        
        Encapsulates capabilities needed from both Epetra and mpi4py
        """
        
        def __init__(self, epetra_comm, mpi4py_comm):
            self.epetra_comm = epetra_comm
            self.mpi4py_comm = mpi4py_comm
            
        @property
        def procID(self):
            return self.epetra_comm.MyPID()
            
        @property
        def Nproc(self):
            return self.epetra_comm.NumProc()
            
        def all(self, a, axis=None):
            return self.mpi4py_comm.allreduce(a.all(axis=axis), op=MPI.LAND)

        def any(self, a, axis=None):
            return self.mpi4py_comm.allreduce(a.any(axis=axis), op=MPI.LOR)

        def allclose(self, a, b, rtol=1.e-5, atol=1.e-8):
            return self.mpi4py_comm.allreduce(numerix.allclose(a, b, rtol=rtol, atol=atol), op=MPI.LAND)
            
        def allequal(self, a, b):
            return self.mpi4py_comm.allreduce(numerix.allequal(a, b), op=MPI.LAND)
            
        def bcast(self, obj=None, root=0):
            return self.mpi4py_comm.bcast(obj=obj, root=root)
            
        def allgather(self, sendobj=None, recvobj=None):
            return self.mpi4py_comm.allgather(sendobj=sendobj, recvobj=recvobj)
            
        def sum(self, a, axis=None):
            return self.epetra_comm.SumAll(a.sum(axis=axis))


    parallel = CommWrapper(epetra_comm=Epetra.PyComm(), mpi4py_comm=MPI.COMM_WORLD)
    serial = CommWrapper(epetra_comm=Epetra.SerialComm(), mpi4py_comm=MPI.COMM_SELF)

except ImportError:
    class DummyComm(object):
        @property
        def procID(self):
            return 0
            
        @property
        def Nproc(self):
            return 1
            
        def Barrier(self):
            pass
            
        def all(self, a, axis=None):
            return a.all(axis=axis)

        def any(self, a, axis=None):
            return a.any(axis=axis)

        def allclose(self, a, b, rtol=1.e-5, atol=1.e-8):
            return numerix.allclose(a, b, rtol=rtol, atol=atol)
            
        def allequal(self, a, b):
            return numerix.allequal(a, b)
            
        def bcast(self, obj=None, root=0):
            return obj
            
        def allgather(self, sendobj=None, recvobj=None):
            if recvobj is not None:
                recvobj[:] = sendobj
            else:
                recvobj = sendobj
                
            return recvobj
            
        def sum(self, a, axis=None):
            return a.sum(axis=axis)


            
    parallel = DummyComm()
    serial = DummyComm()

import dump
import numerix
import vector
from dimensions.physicalField import PhysicalField
from numerix import *
from vitals import Vitals

