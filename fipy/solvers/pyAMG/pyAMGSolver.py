#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "pyAMGSolver.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
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

from fipy.solvers.solver import Solver
import os
from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
from fipy.tools import numerix
        
class PyAMGSolver(Solver):
    """
    The base `PyAMGSolver` class.
    
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    @property
    def _matrixClass(self):
        return _ScipyMeshMatrix

    def _solve_(self, L, x, b):
        """
        Establishes a `pyamg.multilevel.multilevel_solver` object based on
        `self.solveFnc` and then solves, populating `relres` with
        a list of residuals.

        :Parameters:
            - `L`: A `fipy.matrices.pysparseMatrix._PysparseMeshMatrix`.
            - `x`: A `numpy.ndarray`.
            - `b`: A `numpy.ndarray`.
        """

        A = L.asformat('csr')
        if self.preconditioner is None:
            M = None
        else:
            M = self.preconditioner._applyToMatrix(A)

        x, info = self.solveFnc(A, b, x, tol=self.tolerance, maxiter=self.iterations, M=M)

        if os.environ.has_key('FIPY_VERBOSE_SOLVER'):
            from fipy.tools.debug import PRINT
            PRINT('info:', info)

        return x
                        
    def _solve(self):

        if self.var.mesh.communicator.Nproc > 1:
            raise Exception("PyAMG solvers cannot be used with multiple processors")
        
        self.var[:] = self._solve_(self.matrix, self.var.value, numerix.array(self.RHSvector))

