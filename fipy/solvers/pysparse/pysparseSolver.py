#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "pysparseSolver.py"
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
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools.pysparseMatrix import _PysparseMatrix
from fipy.solvers.solver import Solver
from pysparse import precon

class PysparseSolver(Solver):
    """
    The base `pysparseSolver` class.
    
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    def __init__(self, *args, **kwargs):
        if self.__class__ is PysparseSolver:
            raise NotImplementedError, "can't instantiate abstract base class"
            
        Solver.__init__(self, *args, **kwargs)

    def _getMatrixClass(self):
        return _PysparseMatrix
    
        
    def _setupPreconditioner(self, A, defaultPrecon=precon.ssor):
        """
        Internal function for setting up a PySparse preconditioner.
        Returns the preconditioning matrix, `P`, and the resulting
        `A` matrix.
        """
        if self.preconditioner is None:
            self.preconditioner = defaultPrecon

        if self.preconditioner.__name__ == "ssor":
            A = A.to_sss()
            P = self.preconditioner(A)
        else:
            P = self.preconditioner(A)
            A = A.to_csr()

        return P, A
         
