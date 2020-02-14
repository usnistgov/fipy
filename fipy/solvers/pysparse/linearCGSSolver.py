from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

import sys

from pysparse.itsolvers import krylov

from fipy.solvers.pysparse.pysparseSolver import PysparseSolver

__all__ = ["LinearCGSSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

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
        Parameters
        ----------
        precon : ~fipy.solvers.pysparse.preconditioners.preconditioner.Preconditioner, optional
        """
        import warnings
        warnings.warn("The Pysparse CGS solver may return incorrect results for some matrices", UserWarning)
        super(LinearCGSSolver, self).__init__(precon=precon, *args, **kwargs)
        self.solveFnc = krylov.cgs
