from __future__ import division
from __future__ import unicode_literals
from builtins import range
__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO
from PyTrilinos import Epetra
from PyTrilinos import Amesos

from fipy.solvers.trilinos.trilinosSolver import TrilinosSolver

__all__ = ["LinearLUSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearLUSolver(TrilinosSolver):

    """
    The `LinearLUSolver` is an interface to the Amesos KLU solver in Trilinos.

    """

    def __init__(self, tolerance=1e-10, criterion="default", precon=None,
                 iterations=10):
        """
        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'unscaled', 'RHS', 'matrix', 'initial'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon
            *ignored*
        """

        TrilinosSolver.__init__(self, tolerance=tolerance, criterion=criterion,
                                iterations=iterations, precon=None)

        if precon is not None:
            import warnings
            warnings.warn("Trilinos KLU solver does not accept preconditioners.",
                           UserWarning, stacklevel=2)
        self.Factory = Amesos.Factory()

    def _adaptDefaultTolerance(self, L, x, b):
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

        tolerance_factor, _ = self._adaptTolerance(L, x, b)

        for iteration in range(self.iterations):
            residualVector, residual = self._residualVectorAndNorm(L, x, b)

            if residual <= self.tolerance * tolerance_factor:
                break

            xError = Epetra.Vector(L.RowMap())

            Problem = Epetra.LinearProblem(L, xError, residualVector)
            Solver = self.Factory.Create(text_to_native_str("Klu"), Problem)
            Solver.Solve()

            x[:] = x - xError

        self._setConvergence(suite="trilinos",
                             code=AztecOO.AZ_normal,
                             iterations=iteration+1,
                             residual=float(residual))

        self.convergence.warn()

