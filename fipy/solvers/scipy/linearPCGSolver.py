from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from scipy.sparse.linalg import cg

from fipy.solvers.scipy.scipyKrylovSolver import _ScipyKrylovSolver

__all__ = ["LinearPCGSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearPCGSolver(_ScipyKrylovSolver):
    """
    The `LinearPCGSolver` is an interface to the CG solver in Scipy,
    with no preconditioning by default.
    """

    solveFnc = staticmethod(cg)

    def _canSolveAsymmetric(self):
        return False
