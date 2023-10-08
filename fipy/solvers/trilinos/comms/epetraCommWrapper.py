#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "epetraCommWrapper.py"
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

from PyTrilinos import Epetra

from fipy.tools import numerix
from fipy.tools.comms.commWrapper import CommWrapper

__all__ = ["EpetraCommWrapper"]

class EpetraCommWrapper(CommWrapper):
    """MPI Communicator wrapper
    
    Encapsulates capabilities needed for Epetra. 
    Some capabilities are not parallel.
    """
    
    def __init__(self):
        self.epetra_comm = Epetra.PyComm()
        super(EpetraCommWrapper, self).__init__() 
        
    @property
    def procID(self):
        return self.epetra_comm.MyPID()
        
    @property
    def Nproc(self):
        return self.epetra_comm.NumProc()
        
    def Barrier(self):
        self.epetra_comm.Barrier()

    def sum(self, a, axis=None):
        summed = numerix.array(a).sum(axis=axis)
        shape = summed.shape
        if shape == ():
            summed = summed.reshape((1,))
        parallelSummed = self.epetra_comm.SumAll(summed)
        if shape == ():
            parallelSummed = parallelSummed.reshape(())
        return parallelSummed

    def __setstate__(self, dict):
        self.__init__()
        
    def Norm2(self, vec):
        return vec.Norm2()
        
    def MaxAll(self, vec):
        return self.epetra_comm.MaxAll(numerix.array(vec))
        
    def MinAll(self, vec):
        return self.epetra_comm.MinAll(numerix.array(vec))
