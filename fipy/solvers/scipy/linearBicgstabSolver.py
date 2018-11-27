


__docformat__ = 'restructuredtext'

from fipy.solvers.scipy.scipyKrylovSolver import _ScipyKrylovSolver
from scipy.sparse.linalg import bicgstab

__all__ = ["LinearBicgstabSolver"]

class LinearBicgstabSolver(_ScipyKrylovSolver):
    """
    The `LinearBicgstabSolver` is an interface to the Bicgstab solver in
    Scipy, with no preconditioning by default.
    """

    def __init__(self, tolerance=1e-15, iterations=2000, precon=None):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner to use.
        """

        super(LinearBicgstabSolver, self).__init__(tolerance=tolerance, iterations=iterations, precon=precon)
        self.solveFnc = bicgstab
