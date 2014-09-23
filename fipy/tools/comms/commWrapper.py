#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "commWrapper.py"
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
     
    def bcast(self, obj=None, root=0):
        return obj
     
    def allgather(self, sendobj=None, recvobj=None):
        if recvobj is not None:
            recvobj[:] = sendobj
        else:
            recvobj = sendobj
         
        return recvobj
                    
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
