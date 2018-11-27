


__docformat__ = 'restructuredtext'

from scipy.sparse.linalg import gmres

from fipy.solvers.scipy.scipyKrylovSolver import _ScipyKrylovSolver

__all__ = ["LinearGMRESSolver"]

class LinearGMRESSolver(_ScipyKrylovSolver):
    """
    The `LinearGMRESSolver` is an interface to the GMRES solver in
    Scipy, with no preconditioning by default.
    """

    def __init__(self, tolerance=1e-15, iterations=2000, precon=None):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner to use.
        """

        super(LinearGMRESSolver, self).__init__(tolerance=tolerance, iterations=iterations, precon=precon)
        self.solveFnc = gmres
