from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from builtins import range
__docformat__ = 'restructuredtext'

from fipy.solvers.pysparse.pysparseSolver import PysparseSolver
from fipy.matrices.pysparseMatrix import _PysparseMatrixFromShape

__all__ = ["LinearJORSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearJORSolver(PysparseSolver):
    """

    The `LinearJORSolver` solves a linear system of equations using
    Jacobi over-relaxation. This method solves systems with a general
    non-symmetric coefficient matrix.

    """

    def __init__(self, tolerance="default", criterion="default",
                 iterations="default", relaxation=1.0):
        """
        Create a `LinearJORSolver` object.

        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'unscaled', 'RHS', 'matrix', 'initial', 'legacy'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        relaxation : float
            Fraction of update to apply
        """
        super(LinearJORSolver, self).__init__(tolerance=tolerance, criterion=criterion,
                                              iterations=iterations, precon=None)
        self.relaxation = relaxation

    def _solve_(self, L, x, b):

        d = L.takeDiagonal()
        D = _PysparseMatrixFromShape(size=len(d))
        D.putDiagonal(d)

        LU = L - D
        tol = 1e+10
        xold = x.copy()

        tolerance_scale, _ = self._adaptTolerance(L, x, b)

        for iteration in range(self.iterations):
            residual = numerix.L2norm(L * x - b)

            if residual <= self.tolerance * tolerance_scale:
                break

            xold[:] = x
            x[:] = (-(LU) * x + b) / d

            x[:] = xold + self.relaxation * (x - xold)

        self._setConvergence(suite="pysparse",
                             code=0,
                             iterations=iteration+1,
                             residual=residual / tolerance_scale)
