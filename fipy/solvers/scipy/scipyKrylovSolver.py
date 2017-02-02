from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

import os
import warnings

from fipy.solvers.scipy.scipySolver import _ScipySolver
from fipy.solvers import (MaximumIterationWarning,
                          IllegalInputOrBreakdownWarning)

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
                                callback=self._countIterations)
                                
        self.status['iterations'] = self.actualIterations
        if info == 0:
            self.status['code'] = "Success"
        elif info < 0:
            self.status['code'] = IllegalInputOrBreakdownWarning.__class__.__name__
        elif info > 0:
            self.status['code'] = MaximumIterationWarning.__class__.__name__
            
        self._raiseWarning(info, self.actualIterations, 0.)

        if info < 0:
            self._log.debug('failure: %s', self._warningList[info].__class__.__name__)

        return x

    def _raiseWarning(self, info, iter, relres):
        # 0 : successful exit
        # >0 : convergence to tolerance not achieved, number of iterations
        # <0 : illegal input or breakdown

        if info < 0:
            # is stacklevel=5 always what's needed to get to the user's scope?
            warnings.warn(IllegalInputOrBreakdownWarning(self, iter, relres), stacklevel=5)
        elif info > 0:
            warnings.warn(MaximumIterationWarning(self, iter, relres), stacklevel=5)
