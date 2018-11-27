


__docformat__ = 'restructuredtext'

from scipy.sparse.linalg import cg

from fipy.solvers.scipy.scipyKrylovSolver import _ScipyKrylovSolver

__all__ = ["LinearPCGSolver"]

class LinearPCGSolver(_ScipyKrylovSolver):
    """
    The `LinearPCGSolver` is an interface to the CG solver in Scipy,
    with no preconditioning by default.
    """

    def __init__(self, tolerance=1e-15, iterations=2000, precon=None):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner to use.
        """

        super(LinearPCGSolver, self).__init__(tolerance=tolerance, iterations=iterations, precon=precon)
        self.solveFnc = cg

    def _canSolveAsymmetric(self):
        return False
