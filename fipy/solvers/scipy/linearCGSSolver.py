from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from scipy.sparse.linalg import cgs

from fipy.solvers.scipy.scipyKrylovSolver import _ScipyKrylovSolver

__all__ = ["LinearCGSSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearCGSSolver(_ScipyKrylovSolver):
    """
    The `LinearCGSSolver` is an interface to the CGS solver in Scipy,
    with no preconditioning by default.
    """

    solveFnc = staticmethod(cgs)
