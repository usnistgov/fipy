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

    def __init__(self, tolerance=1e-10, iterations=10, precon=None, maxIterations=10):
        """
        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        iterations : int
            Maximum number of iterative steps to perform.
        precon
            *ignored*
        """

        iterations = min(iterations, maxIterations)

        TrilinosSolver.__init__(self, tolerance=tolerance,
                                iterations=iterations, precon=None)

        if precon is not None:
            import warnings
            warnings.warn("Trilinos KLU solver does not accept preconditioners.",
                           UserWarning, stacklevel=2)
        self.Factory = Amesos.Factory()


    def _solve_(self, L, x, b):

        for iteration in range(self.iterations):
             # residualVector = L*x - b
             residualVector = Epetra.Vector(L.RangeMap())
             L.Multiply(False, x, residualVector)
             # If A is an Epetra.Vector with map M
             # and B is an Epetra.Vector with map M
             # and C = A - B
             # then C is an Epetra.Vector with *no map* !!!?!?!
             residualVector -= b

             residual = residualVector.Norm2()
             if iteration == 0:
                 residual0 = residual

             if residual <= self.tolerance * residual0:
                 break

             xError = Epetra.Vector(L.RowMap())

             Problem = Epetra.LinearProblem(L, xError, residualVector)
             Solver = self.Factory.Create(text_to_native_str("Klu"), Problem)
             Solver.Solve()

             x[:] = x - xError

        self._setConvergence(suite="trilinos"
                             code=AztecOO.AZ_normal,
                             iterations=iteration+1,
                             residual=residual)

        self.convergence.warn()

