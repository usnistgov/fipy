#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "parallelEpetraCommWrapper.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

from PyTrilinos import Epetra
from mpi4py import MPI

from fipy.tools import numerix
from fipy.solvers.trilinos.comms.epetraCommWrapper import EpetraCommWrapper

__all__ = ["ParallelEpetraCommWrapper"]

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

    def bcast(self, obj=None, root=0):
        return self.mpi4py_comm.bcast(obj=obj, root=root)

    def allgather(self, sendobj=None, recvobj=None):
        return self.mpi4py_comm.allgather(sendobj=sendobj, recvobj=recvobj)
