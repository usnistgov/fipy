from __future__ import division
from __future__ import unicode_literals
from builtins import range
__docformat__ = 'restructuredtext'

from scipy.sparse.linalg import splu

from fipy.solvers.scipy.scipySolver import _ScipySolver
from fipy.tools import numerix
from fipy.tools.timer import Timer

__all__ = ["LinearLUSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearLUSolver(_ScipySolver):
    """
    The `LinearLUSolver` solves a linear system of equations using
    LU-factorization.  The `LinearLUSolver` is a wrapper class for the
    the Scipy `scipy.sparse.linalg.splu` module.
    """

    def _solve_(self, L, x, b):
        diag = L.takeDiagonal()
        maxdiag = max(numerix.absolute(diag))

        self._log.debug("BEGIN precondition")

        with Timer() as t:
            L = L * (1 / maxdiag)
            b = b * (1 / maxdiag)

        self._log.debug("END precondition - {} ns".format(t.elapsed))

        self._log.debug("BEGIN solve")

        with Timer() as t:
            LU = splu(L.matrix.asformat("csc"), diag_pivot_thresh=1.,
                                                relax=1,
                                                panel_size=10,
                                                permc_spec=3)

            error0 = numerix.sqrt(numerix.sum((L * x - b)**2))

            for iteration in range(min(self.iterations, 10)):
                errorVector = L * x - b

                if numerix.sqrt(numerix.sum(errorVector**2))  <= self.tolerance * error0:
                    break

                xError = LU.solve(errorVector)
                x[:] = x - xError

        self._log.debug("END solve - {} ns".format(t.elapsed))

        self._log.debug('iterations: %d / %d', iteration+1, self.iterations)
        self._log.debug('residual: %s', numerix.sqrt(numerix.sum(errorVector**2)))

        return x
