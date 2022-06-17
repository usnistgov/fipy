from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

import logging

_log = logging.getLogger(__name__)

import os

from fipy.solvers.scipy.scipySolver import _ScipySolver

class _ScipyKrylovSolver(_ScipySolver):
    """
    The base `ScipyKrylovSolver` class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def _solve_(self, L, x, b):
        A = L.matrix
        if self.preconditioner is None:
            M = None
        else:
            M = self.preconditioner._applyToMatrix(A)

        x, info = self.solveFnc(A, b, x,
                                tol=self.tolerance,
                                maxiter=self.iterations,
                                M=M,
                                atol='legacy')

        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            from fipy.tools.debug import PRINT
            if info < 0:
                PRINT('failure', self._warningList[info].__class__.__name__)

        if info < 0:
            _log.debug('%s failure: %s', self._warningList[info].__class__.__name__)

        return x
