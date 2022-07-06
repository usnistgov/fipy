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

    def _solve_(self, L, x, b):
        A = L.matrix
        if self.preconditioner is None:
            M = None
        else:
            M = self.preconditioner._applyToMatrix(A)

        self.actualIterations = 0
        x, info = self.solveFnc(A, b, x,
                                tol=self.tolerance,
                                maxiter=self.iterations,
                                M=M,
                                atol='legacy',
                                callback=self._countIterations,
                                callback_type='legacy')

        self._setConvergence(suite="scipy",
                             code=numerix.sign(info),
                             actual_code=info,
                             iterations=self.actualIterations,
                             residual=0.)

        self.convergence.warn()

        return x
