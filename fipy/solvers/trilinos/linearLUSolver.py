from __future__ import division
from builtins import range
from past.utils import old_div
__docformat__ = 'restructuredtext'

import os

from PyTrilinos import Epetra
from PyTrilinos import Amesos

from fipy.solvers.trilinos.trilinosSolver import TrilinosSolver

__all__ = ["LinearLUSolver"]

class LinearLUSolver(TrilinosSolver):

    """
    The `LinearLUSolver` is an interface to the Amesos KLU solver in Trilinos.

    """

    def __init__(self, tolerance=1e-10, iterations=10, precon=None, maxIterations=10):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.

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
             # errorVector = L*x - b
             errorVector = Epetra.Vector(L.RangeMap())
             L.Multiply(False, x, errorVector)
             # If A is an Epetra.Vector with map M
             # and B is an Epetra.Vector with map M
             # and C = A - B
             # then C is an Epetra.Vector with *no map* !!!?!?!
             errorVector -= b

             tol = errorVector.Norm1()

             if iteration == 0:
                 tol0 = tol

             if (old_div(tol, tol0)) <= self.tolerance:
                 break

             xError = Epetra.Vector(L.RowMap())

             Problem = Epetra.LinearProblem(L, xError, errorVector)
             Solver = self.Factory.Create("Klu", Problem)
             Solver.Solve()

             x[:] = x - xError

        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            from fipy.tools.debug import PRINT
            PRINT('iterations: %d / %d' % (iteration + 1, self.iterations))
            PRINT('residual:', errorVector.Norm2())
