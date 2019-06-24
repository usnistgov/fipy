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

    def __init__(self, tolerance=1e-15, iterations=2000, precon=None):
        """
        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        iterations : int
            Maximum number of iterative steps to perform.
        precon
            Preconditioner to use.
        """

        super(LinearGMRESSolver, self).__init__(tolerance=tolerance, iterations=iterations, precon=precon)
        self.solveFnc = gmres
