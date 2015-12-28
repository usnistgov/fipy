from mpi4py import MPI

m4comm = MPI.COMM_WORLD
mpi4py_info = "mpi4py: processor %d of %d" % (m4comm.Get_rank(),
                                              m4comm.Get_size())

from PyTrilinos import Epetra

epcomm = Epetra.PyComm()

trilinos_info = "PyTrilinos: processor %d of %d" % (epcomm.MyPID(),
                                                    epcomm.NumProc())

from fipy import parallelComm, Grid1D

mesh = Grid1D(nx=10)

fipy_info = "FiPy: %d cells on processor %d of %d" % (mesh.numberOfCells,
                                                      parallelComm.procID,
                                                      parallelComm.Nproc)

print " :: ".join((mpi4py_info, trilinos_info, fipy_info))
