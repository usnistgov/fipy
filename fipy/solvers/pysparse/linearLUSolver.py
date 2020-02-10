from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from builtins import range
__docformat__ = 'restructuredtext'

import os

from pysparse.direct import superlu

from fipy.solvers.pysparse.pysparseSolver import PysparseSolver
from fipy.tools import numerix

DEBUG = False

__all__ = ["LinearLUSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

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

        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.pysparse.preconditioners.preconditioner.Preconditioner
            *ignored*
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
