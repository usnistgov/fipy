__docformat__ = 'restructuredtext'

import sys

from pysparse import itsolvers

from fipy.solvers.pysparse.pysparseSolver import PysparseSolver

__all__ = ["LinearCGSSolver"]

class LinearCGSSolver(PysparseSolver):

    """

    The `LinearCGSSolver` solves a linear system of equations using
    the conjugate gradient squared method (CGS), a variant of the
    biconjugate gradient method (BiCG). CGS solves linear systems with
    a general non-symmetric coefficient matrix.

    The `LinearCGSSolver` is a wrapper class for the the Pysparse_
    `itsolvers.cgs()` method.

    .. _Pysparse: http://pysparse.sourceforge.net

    """

    def __init__(self, precon=None, *args, **kwargs):
        """
        :Parameters:
          - `precon`: Preconditioner to use
        """
        import warnings
        warnings.warn("The Pysparse CGS solver may return incorrect results for some matrices", UserWarning)
        super(LinearCGSSolver, self).__init__(precon=precon, *args, **kwargs)
        self.solveFnc = itsolvers.cgs
