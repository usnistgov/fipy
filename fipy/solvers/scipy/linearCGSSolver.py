


__docformat__ = 'restructuredtext'

from scipy.sparse.linalg import cgs

from fipy.solvers.scipy.scipyKrylovSolver import _ScipyKrylovSolver

__all__ = ["LinearCGSSolver"]

class LinearCGSSolver(_ScipyKrylovSolver):
    """
    The `LinearCGSSolver` is an interface to the CGS solver in Scipy,
    with no preconditioning by default.
    """

    def __init__(self, tolerance=1e-15, iterations=2000, precon=None):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner to use.
        """

        super(LinearCGSSolver, self).__init__(tolerance=tolerance, iterations=iterations, precon=precon)
        self.solveFnc = cgs
