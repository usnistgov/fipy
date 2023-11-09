from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from pysparse.itsolvers import krylov

from .linearRHSSolver import LinearRHSSolver
from .preconditioners import SSORPreconditioner

__all__ = ["LinearPCGSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearPCGSolver(LinearRHSSolver):
    """

    The `LinearPCGSolver` solves a linear system of equations using the
    preconditioned conjugate gradient method (PCG) with symmetric successive
    over-relaxation (SSOR) preconditioning by default. Alternatively,
    Jacobi preconditioning can be specified through `precon`.
    The PCG method solves systems with
    a symmetric positive definite coefficient matrix.

    The `LinearPCGSolver` is a wrapper class for the the Pysparse_
    `itsolvers.pcg()` and `precon.ssor()` methods.

    .. _Pysparse: http://pysparse.sourceforge.net

    """

    solveFnc = staticmethod(krylov.pcg)

    DEFAULT_PRECONDITIONER = SSORPreconditioner

    def _canSolveAsymmetric(self):
        return False
