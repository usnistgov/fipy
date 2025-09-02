from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from builtins import range
__docformat__ = 'restructuredtext'

from pysparse.direct import superlu

from .pysparseSolver import PysparseSolver
from fipy.tools import numerix
from fipy.tools.timer import Timer

__all__ = ["LinearLUSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearLUSolver(PysparseSolver):
    """Interface to :term:`LU`-factorization in :ref:`Pysparse`.
    """

    def __init__(self, tolerance="default", criterion="default",
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
        precon
            *ignored*
        """
        super(LinearLUSolver, self).__init__(tolerance=tolerance, criterion=criterion,
                                             iterations=iterations, precon=None)

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
        """Solve system of equations posed for PySparse

        Parameters
        ----------
        L : ~pysparse.spmatrix.ll_mat
            Sparse matrix
        x : array_like
            Solution vector
        b : array_like
            Right hand side vector

        Returns
        -------
        x : ndarray
            Solution vector
        """
        tolerance_scale, _ = self._adaptTolerance(L, x, b)

        self._log.debug("BEGIN solve")

        with Timer() as t:
            LU = superlu.factorize(L.to_csr())

            for iteration in range(self.iterations):
                residualVector, residual = self._residualVectorAndNorm(L, x, b)

                if residual <= self.tolerance * tolerance_scale:
                    break

                xError = numerix.zeros(len(b), 'd')

                LU.solve(residualVector, xError)
                x[:] = x - xError

        self._log.debug("END solve - {} ns".format(t.elapsed))

        self._setConvergence(suite="pysparse",
                             code=0,
                             iterations=iteration+1,
                             residual=residual)

        self.convergence.warn()

        return x
