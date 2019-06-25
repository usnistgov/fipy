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
    def __init__(self, tolerance=1e-10, iterations=1000, relaxation=1.0):
        """
        The `Solver` class should not be invoked directly.

        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        iterations : int
            Maximum number of iterative steps to perform.
        relaxation : float
            Fraction of update to apply
        """
        super(LinearJORSolver, self).__init__(tolerance=tolerance,
                                              iterations=iterations)
        self.relaxation = relaxation

    def _solve_(self, L, x, b):

        d = L.takeDiagonal()
        D = _PysparseMatrixFromShape(size=len(d))
        D.putDiagonal(d)

        LU = L - D
        tol = 1e+10
        xold = x.copy()

        for iteration in range(self.iterations):
            if tol <= self.tolerance:
                break

            residual = L * x - b

            xold[:] = x
            x[:] = (-(LU) * x + b) / d

            x[:] = xold + self.relaxation * (x - xold)

            tol = max(abs(residual))

            print(iteration, tol)
