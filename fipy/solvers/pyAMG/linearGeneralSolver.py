


from fipy.solvers.scipy.scipySolver import _ScipySolver
from pyamg import solve
import os
from fipy.tools import numerix

__all__ = ["LinearGeneralSolver"]

class LinearGeneralSolver(_ScipySolver):
    """
    The `LinearGeneralSolver` is an interface to the generic pyAMG,
    which solves the arbitrary system Ax=b with the best out-of-the box
    choice for a solver. See `pyAMG.solve` for details.
    """

    def _solve_(self, L, x, b):
        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            verbosity = True
        else:
            verbosity = False

        return solve(L.matrix, b, verb=verbosity, tol=self.tolerance)
