#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "linearGeneralSolver.py"
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

from fipy.solvers.pyAMG.pyAMGSolver import PyAMGSolver
from scipy.sparse.linalg import gmres
from pyamg import solve
import os

class LinearGeneralSolver(PyAMGSolver):

    def __init__(self, tolerance=1e-15, iterations=2000, steps=None, precon=None):
        PyAMGSolver.__init__(self, tolerance=tolerance, iterations=iterations, steps=steps, precon=precon)
        
    def _canSolveAssymetric(self):
        return True

    def _solve_(self, L, x, b):
        """

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

        if os.environ.has_key('FIPY_VERBOSE_SOLVER'):
            verbosity = True
        else:
            verbosity = False

        return solve(A, b, verb=verbosity, tol=self.tolerance)
