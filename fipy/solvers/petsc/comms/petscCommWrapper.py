#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "petscCommWrapper.py"
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

from petsc4py import PETSc

from fipy.tools import numerix
from fipy.tools.comms.commWrapper import CommWrapper

__all__ = ["PETScCommWrapper"]

class PETScCommWrapper(CommWrapper):
    """MPI Communicator wrapper
    
    Encapsulates capabilities needed for PETSc. 
    Some capabilities are not parallel.
    """
    
    def __init__(self, petsc4py_comm=PETSc.COMM_WORLD):
        self.petsc4py_comm = petsc4py_comm
        super(PETScCommWrapper, self).__init__()
        
    @property
    def mpi4py_comm(self):
        return self.petsc4py_comm.tompi4py()
        
    def Norm2(self, vec):
        return vec.norm(norm_type=1)
