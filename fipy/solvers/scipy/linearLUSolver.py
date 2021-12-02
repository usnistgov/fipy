from __future__ import division
from __future__ import unicode_literals
from builtins import range
__docformat__ = 'restructuredtext'

import os

from scipy.sparse.linalg import splu

from fipy.solvers.scipy.scipySolver import _ScipySolver
from fipy.tools import numerix

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

        L = L * (1 / maxdiag)
        b = b * (1 / maxdiag)

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

        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            from fipy.tools.debug import PRINT
            PRINT('iterations: %d / %d' % (iteration+1, self.iterations))
            PRINT('residual:', numerix.sqrt(numerix.sum(errorVector**2)))

        return x
