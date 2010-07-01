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
            raise NotImplementedError, \
                  "can't instantiate abstract base class"
            
        Solver.__init__(self, *args, **kwargs)

    def _getMatrixClass(self):
        return _PysparseMatrix

    def _solve_(self, L, x, b):
        """
        `_solve_` is only for use by solvers which may use
        preconditioning. If you are writing a solver which
        doesn't use preconditioning, this must be overridden.
        """
        A = L._getMatrix()

        if self.preconditioner is None:
            P = None
        else:
            P, A = self.preconditioner._applyToMatrix(A)

        info, iter, relres = self.solveFnc(A, b, x, self.tolerance, 
                                           self.iterations, P)
        
        self._raiseWarning(info, iter, relres)
          
    def _solve(self):

        if self.var.getMesh().parallelModule.Nproc > 1:
            raise Exception("PySparse solvers cannot be used with multiple processors")
        
        array = self.var.getNumericValue()

        self._solve_(self.matrix, array, self.RHSvector)
        factor = self.var.getUnit().factor
        if factor != 1:
            array /= self.var.getUnit().factor
        self.var[:] = array
