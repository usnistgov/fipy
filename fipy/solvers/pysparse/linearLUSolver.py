from __future__ import print_function
__docformat__ = 'restructuredtext'

import os

from pysparse import superlu

from fipy.solvers.pysparse.pysparseSolver import PysparseSolver
from fipy.tools import numerix

DEBUG = False

__all__ = ["LinearLUSolver"]

class LinearLUSolver(PysparseSolver):
    """

    The `LinearLUSolver` solves a linear system of equations using
    LU-factorization. This method solves systems with a general
    non-symmetric coefficient matrix using partial pivoting.

    The `LinearLUSolver` is a wrapper class for the the Pysparse_
    `superlu.factorize()` method.

    .. _Pysparse: http://pysparse.sourceforge.net

    """

    def __init__(self, tolerance=1e-10, iterations=10,
                       maxIterations=10, precon=None):
        """
        Creates a `LinearLUSolver`.

        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The number of LU decompositions to perform.
            For large systems a number of iterations is generally required.
          - `precon`: not used but maintains a common interface.

        """

        iterations = min(iterations, maxIterations)

        super(LinearLUSolver, self).__init__(tolerance = tolerance,
                                             iterations = iterations)

    def _solve_(self, L, x, b):
        diag = L.takeDiagonal()
        maxdiag = max(numerix.absolute(diag))

        L = L * (1 / maxdiag)
        b = b * (1 / maxdiag)

        LU = superlu.factorize(L.matrix.to_csr())

        if DEBUG:
            import sys
            print(L.matrix, file=sys.stderr)

        error0 = numerix.sqrt(numerix.sum((L * x - b)**2))

        for iteration in range(self.iterations):
            errorVector = L * x - b

            if (numerix.sqrt(numerix.sum(errorVector**2)) / error0)  <= self.tolerance:
                break

            xError = numerix.zeros(len(b), 'd')
            LU.solve(errorVector, xError)
            x[:] = x - xError

        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            from fipy.tools.debug import PRINT
            PRINT('iterations: %d / %d' % (iteration+1, self.iterations))
            PRINT('residual:', numerix.sqrt(numerix.sum(errorVector**2)))
