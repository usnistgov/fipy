from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from ..pysparseMatrixSolver import PysparseMatrixSolver
from fipy.tools import numerix

__all__ = ["PysparseSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class PysparseSolver(PysparseMatrixSolver):
    """
    The base `pysparseSolver` class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

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

        # Pysparse returns the relative residual,
        # which changes depending on which solver is used
        legacy_norm = self._legacyNorm(L, x, b)

        self._log.debug("BEGIN precondition")

        if self.preconditioner is None:
            P = None
        else:
            P, L = self.preconditioner._applyToMatrix(L)

        self._log.debug("END precondition")
        self._log.debug("BEGIN solve")

        info, iter, relres = self.solveFnc(L, b, x,
                                           self.tolerance * tolerance_scale,
                                           self.iterations, P)

        self._log.debug("END solve")

        self._setConvergence(suite="pysparse",
                             code=info,
                             iterations=iter,
                             tolerance_scale=tolerance_scale,
                             residual=relres * legacy_norm)

        self.convergence.warn()

        return x

    def _rhsNorm(self, L, x, b):
        return numerix.L2norm(b)

    def _matrixNorm(self, L, x, b):
        return L.norm('inf')

    def _residualVectorAndNorm(self, L, x, b):
        y = numerix.empty((L.shape[0],))
        L.matvec(x, y)
        residualVector = y - b

        return residualVector, numerix.L2norm(residualVector)

    def _adaptUnscaledTolerance(self, L, x, b):
        factor = 1. / self._legacyNorm(L, x, b)
        return (factor, None)

    def _adaptRHSTolerance(self, L, x, b):
        factor = self._rhsNorm(L, x, b) / self._legacyNorm(L, x, b)
        return (factor, None)

    def _adaptMatrixTolerance(self, L, x, b):
        factor = self._matrixNorm(L, x, b) / self._legacyNorm(L, x, b)
        return (factor, None)

    def _adaptInitialTolerance(self, L, x, b):
        factor = self._residualNorm(L, x, b) / self._legacyNorm(L, x, b)
        return (factor, None)

    def _adaptLegacyTolerance(self, L, x, b):
        return (1., None)
