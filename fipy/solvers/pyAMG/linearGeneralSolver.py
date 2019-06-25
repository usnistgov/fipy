from __future__ import unicode_literals
from fipy.solvers.scipy.scipySolver import _ScipySolver
from pyamg import solve
import os
from fipy.tools import numerix

__all__ = ["LinearGeneralSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearGeneralSolver(_ScipySolver):
    """
    The `LinearGeneralSolver` is an interface to the generic PyAMG,
    which solves the arbitrary system Ax=b with the best out-of-the box
    choice for a solver. See `pyAMG.solve` for details.
    """

    def _solve_(self, L, x, b):
        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            verbosity = True
        else:
            verbosity = False

        return solve(L.matrix, b, verb=verbosity, tol=self.tolerance)
