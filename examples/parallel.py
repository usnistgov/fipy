from __future__ import print_function
from __future__ import unicode_literals

from fipy import parallelComm

def main():
    titles = []
    results = []

    titles.append("mpi4py")

    try:
        from mpi4py import MPI

        m4comm = MPI.COMM_WORLD
        results.append("processor %d of %d" % (m4comm.Get_rank(),
                                               m4comm.Get_size()))

    except Exception as e:
        results.append(str(e))



    titles.append("PyTrilinos")

    try:
        from PyTrilinos import Epetra

        epcomm = Epetra.PyComm()

        results.append("processor %d of %d" % (epcomm.MyPID(),
                                               epcomm.NumProc()))
    except Exception as e:
        results.append(str(e))


    titles.append("petsc4py")

    try:
        from petsc4py import PETSc

        pecomm = PETSc.COMM_WORLD

        results.append("processor %d of %d" % (pecomm.rank,
                                               pecomm.size))
    except Exception as e:
        results.append(str(e))


    titles.append("FiPy")

    try:
        from fipy import Grid1D

        mesh = Grid1D(nx=10)

        results.append("%d cells on processor %d of %d" % (mesh.numberOfCells,
                                                           parallelComm.procID,
                                                           parallelComm.Nproc))
    except Exception as e:
        results.append(str(e))



    lengths = [parallelComm.MaxAll(len(s)) for s in results]
    formats = ["{{{0}:^{1}}}".format(i, l) for i, l in enumerate(lengths)]

    if parallelComm.procID == 0:
        print("    ".join(formats).format(*titles))

    parallelComm.Barrier()

    print(" :: ".join(formats).format(*results))

if __name__ == "__main__":
    main()
