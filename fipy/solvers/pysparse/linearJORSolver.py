#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "linearJORSolver.py"
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

from fipy.solvers.pysparse.pysparseSolver import PysparseSolver
from fipy.matrices.pysparseMatrix import _PysparseMatrixFromShape

__all__ = ["LinearJORSolver"]

class LinearJORSolver(PysparseSolver):
    """

    The `LinearJORSolver` solves a linear system of equations using
    Jacobi over-relaxation. This method solves systems with a general
    non-symmetric coefficient matrix.

    """
    def __init__(self, tolerance=1e-10, iterations=1000, relaxation=1.0):
        """
        The `Solver` class should not be invoked directly.

        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `relaxation`: The relaxation.

        """
        super(LinearJORSolver, self).__init__(tolerance=tolerance,
                                              iterations=iterations)
        self.relaxation = relaxation

    def _solve_(self, L, x, b):

        d = L.takeDiagonal()
        D = _PysparseMatrixFromShape(size=len(d))
        D.putDiagonal(d)

        LU = L - D
        tol = 1e+10
        xold = x.copy()

        for iteration in range(self.iterations):
            if tol <= self.tolerance:
                break

            residual = L * x - b

            xold[:] = x
            x[:] = (-(LU) * x + b) / d

            x[:] = xold + self.relaxation * (x - xold)

            tol = max(abs(residual))

            print iteration,tol
