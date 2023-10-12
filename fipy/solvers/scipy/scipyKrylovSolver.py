from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

import os
import warnings

from .scipySolver import _ScipySolver
from fipy.tools import numerix

class _ScipyKrylovSolver(_ScipySolver):
    """
    The base `ScipyKrylovSolver` class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def _countIterations(self, xk):
        self.actualIterations += 1

    def _adaptDefaultTolerance(self, L, x, b):
        return self._adaptRHSTolerance(L, x, b)

    def _adaptUnscaledTolerance(self, L, x, b):
        factor = 1. / self._rhsNorm(L, x, b)
        return (factor, None)

    def _adaptRHSTolerance(self, L, x, b):
        return (1., None)

    def _adaptMatrixTolerance(self, L, x, b):
        factor = self._matrixNorm(L, x, b) / self._rhsNorm(L, x, b)
        return (factor, None)

    def _adaptInitialTolerance(self, L, x, b):
        factor = self._residualNorm(L, x, b) / self._rhsNorm(L, x, b)
        return (factor, None)

    def _solve_(self, L, x, b):
        A = L.matrix
        if self.preconditioner is None:
            M = None
        else:
            M = self.preconditioner._applyToMatrix(A)

        tolerance_scale, _ = self._adaptTolerance(L, x, b)

        self.actualIterations = 0

        rtol, atol = (self.scale_tolerance(tol, tolerance_scale)
                      for tol in (self.tolerance,
                                  self.absolute_tolerance))

        self._log.debug("BEGIN solve")

        x, info = self.solveFnc(A, b, x,
                                tol=rtol,
                                atol=atol,
                                maxiter=self.iterations,
                                M=M,
                                callback=self._countIterations)

        self._log.debug("END solve")

        self._setConvergence(suite="scipy",
                             code=numerix.sign(info),
                             actual_code=info,
                             iterations=self.actualIterations,
                             residual=self._residualNorm(L, x, b))

        self.convergence.warn()

        return x
