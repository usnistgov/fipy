#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "mpi4pyCommWrapper.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

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
