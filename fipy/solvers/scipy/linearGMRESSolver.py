from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from scipy.sparse.linalg import gmres

from fipy.solvers.scipy.scipyKrylovSolver import _ScipyKrylovSolver

__all__ = ["LinearGMRESSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearGMRESSolver(_ScipyKrylovSolver):
    """
    The `LinearGMRESSolver` is an interface to the GMRES solver in
    Scipy, with no preconditioning by default.
    """

    solveFnc = staticmethod(gmres)
