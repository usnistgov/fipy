from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from .pysparseSolver import PysparseSolver

class LinearRHSSolver(PysparseSolver):
    """Wrapper for solvers that normalize the residual by the right-hand-side.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def __init__(self, tolerance=1e-10, criterion="legacy",
                 iterations=1000, precon=None):
        """
        Create a `Solver` object.

        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'unscaled', 'RHS', 'matrix', 'initial', 'legacy'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.pysparse.preconditioners.preconditioner.Preconditioner
            Preconditioner to use.
        """
        super(LinearRHSSolver, self).__init__(tolerance=tolerance, criterion=criterion,
                                              iterations=iterations, precon=precon)

    def _legacyNorm(self, L, x, b):
        return self._rhsNorm(L, x, b)

    def _adaptLegacyTolerance(self, L, x, b):
        return self._adaptRHSTolerance(L, x, b)

    def _adaptUnscaledTolerance(self, L, x, b):
        factor = 1. / self._rhsNorm(L, x, b)
        return (factor, None)

    def _adaptRHSTolerance(self, L, x, b):
        return (1., None)

    def _adaptMatrixTolerance(self, L, x, b):
        factor = self._matrixNorm(L, x, b) / self._rhsNorm(L, x, b)
        return (factor, None)

    def _adaptInitialTolerance(self, L, x, b):
        factor = self._residualNorm(L, x, b) / self._rhsNorm(L, x, b)
        return (factor, None)
