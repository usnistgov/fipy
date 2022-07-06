from __future__ import division
from __future__ import unicode_literals
from builtins import range
__docformat__ = 'restructuredtext'

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

        for iteration in range(min(self.iterations, 10)):
            residualVector = L * x - b

            residual = numerix.L2norm(residualVector)
            if iteration == 0:
                residual0 = residual

            if residual <= self.tolerance * residual0:
                break

            xError = LU.solve(residualVector)
            x[:] = x - xError

        self._setConvergence(suite="scipy",
                             code=0,
                             iterations=iteration+1,
                             residual=residual)

        self.convergence.warn()

        return x
