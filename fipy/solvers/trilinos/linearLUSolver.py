from __future__ import division
from __future__ import unicode_literals
from builtins import range
__docformat__ = 'restructuredtext'

from PyTrilinos import Epetra
from PyTrilinos import Amesos

from fipy.solvers.trilinos.trilinosSolver import TrilinosSolver
from fipy.tools.timer import Timer

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

        self._log.debug("BEGIN solve")

        with Timer() as t:
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

                 if (tol / tol0) <= self.tolerance:
                     break

                 xError = Epetra.Vector(L.RowMap())

                 Problem = Epetra.LinearProblem(L, xError, errorVector)
                 Solver = self.Factory.Create(text_to_native_str("Klu"), Problem)
                 Solver.Solve()

                 x[:] = x - xError

        self._log.debug("END solve - {} ns".format(t.elapsed))

        self._log.debug('iterations: %d / %d', iteration+1, self.iterations)
        self._log.debug('residual: %s', errorVector.Norm2())
