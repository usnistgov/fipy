from __future__ import division
from builtins import range
__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO
from PyTrilinos import Epetra
from PyTrilinos import Amesos

from fipy.solvers.trilinos.trilinosSolver import TrilinosSolver
from fipy.tools.timer import Timer

__all__ = ["LinearLUSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearLUSolver(TrilinosSolver):

    """Interface to the Amesos KLU solver in :ref:`TRILINOS`.

    KLU is a direct, serial :term:`LU` solver.
    """

    def __init__(self, tolerance="default", criterion="default", precon=None,
                 iterations=10):
        """
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

        if precon is not None:
            import warnings
            warnings.warn("Trilinos KLU solver does not accept preconditioners.",
                           UserWarning, stacklevel=2)
        self.Factory = Amesos.Factory()

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
        """Solve system of equations posed for PyTrilinos

        Parameters
        ----------
        L : Epetra.CrsMatrix
            Sparse matrix
        x : Epetra.Vector
            Solution variable as non-ghosted vector
        b : Epetra.Vector
            Right-hand side as non-ghosted vector

        Returns
        -------
        x : Epetra.Vector
            Solution variable as non-ghosted vector
        """

        tolerance_scale, _ = self._adaptTolerance(L, x, b)

        self._log.debug("BEGIN solve")

        with Timer() as t:
            for iteration in range(self.iterations):
                residualVector, residual = self._residualVectorAndNorm(L, x, b)

                if residual <= self.tolerance * tolerance_scale:
                    break

                xError = Epetra.Vector(L.RowMap())

                Problem = Epetra.LinearProblem(L, xError, residualVector)
                Solver = self.Factory.Create(text_to_native_str("Klu"), Problem)
                Solver.Solve()

                x[:] = x - xError

        self._log.debug("END solve - {} ns".format(t.elapsed))

        self._setConvergence(suite="trilinos",
                             code=AztecOO.AZ_normal,
                             iterations=iteration+1,
                             residual=float(residual))

        self.convergence.warn()

        return x
