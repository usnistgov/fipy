from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.solvers.scipy.scipySolver import _ScipySolver
from fipy.tools.timer import Timer

class _ScipyKrylovSolver(_ScipySolver):
    """
    The base `ScipyKrylovSolver` class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def _solve_(self, L, x, b):
        A = L.matrix

        self._log.debug("BEGIN precondition")

        with Timer() as t:
            if self.preconditioner is None:
                M = None
            else:
                M = self.preconditioner._applyToMatrix(A)

        self._log.debug("END precondition - {} ns".format(t.elapsed))

        self._log.debug("BEGIN solve")

        with Timer() as t:
            x, info = self.solveFnc(A, b, x,
                                    rtol=self.tolerance,
                                    maxiter=self.iterations,
                                    M=M,
                                    atol=0.0)

        self._log.debug("END solve - {} ns".format(t.elapsed))

        if info < 0:
            self._log.debug('failure: %s', self._warningList[info].__class__.__name__)

        return x
