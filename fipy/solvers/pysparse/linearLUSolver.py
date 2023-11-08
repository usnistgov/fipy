from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from builtins import range
__docformat__ = 'restructuredtext'

from pysparse.direct import superlu

from .pysparseSolver import PysparseSolver
from fipy.tools import numerix

__all__ = ["LinearLUSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearLUSolver(PysparseSolver):
    """

    The `LinearLUSolver` solves a linear system of equations using
    LU-factorization. This method solves systems with a general
    non-symmetric coefficient matrix using partial pivoting.

    The `LinearLUSolver` is a wrapper class for the the Pysparse_
    `superlu.factorize()` method.

    .. _Pysparse: http://pysparse.sourceforge.net

    """

    def __init__(self, tolerance=1e-10, criterion="default",
                 iterations=10, precon=None):
        """
        Creates a `LinearLUSolver`.

        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'unscaled', 'RHS', 'matrix', 'initial', 'legacy'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.pysparse.preconditioners.preconditioner.Preconditioner
            *ignored*
        """
        super(LinearLUSolver, self).__init__(tolerance=tolerance, criterion=criterion,
                                             iterations=iterations)

    def _adaptLegacyTolerance(self, L, x, b):
        return self._adaptInitialTolerance(L, x, b)

    def _adaptUnscaledTolerance(self, L, x, b):
        return (1., None)

    def _adaptRHSTolerance(self, L, x, b):
        return (self._rhsNorm(L, x, b), None)

    def _adaptMatrixTolerance(self, L, x, b):
        return (self._matrixNorm(L, x, b), None)

    def _adaptInitialTolerance(self, L, x, b):
        return (self._residualNorm(L, x, b), None)

    def _solve_(self, L, x, b):
        LU = superlu.factorize(L.matrix.to_csr())

        tolerance_scale, _ = self._adaptTolerance(L, x, b)

        self._log.debug("BEGIN solve")

        for iteration in range(self.iterations):
            residualVector, residual = self._residualVectorAndNorm(L, x, b)

            if residual <= self.tolerance * tolerance_scale:
                break

            xError = numerix.zeros(len(b), 'd')

            LU.solve(residualVector, xError)
            x[:] = x - xError

        self._log.debug("END solve")

        self._setConvergence(suite="pysparse",
                             code=0,
                             iterations=iteration+1,
                             residual=residual)
