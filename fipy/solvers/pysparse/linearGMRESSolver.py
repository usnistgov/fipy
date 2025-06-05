from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from pysparse.itsolvers import krylov

from .linearInitialSolver import LinearInitialSolver
from .preconditioners import JacobiPreconditioner

__all__ = ["LinearGMRESSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearGMRESSolver(LinearInitialSolver):
    """

    The `LinearGMRESSolver` solves a linear system of equations using the
    generalized minimal residual method (GMRES) with Jacobi
    preconditioning. GMRES solves systems with a general non-symmetric
    coefficient matrix.

    The `LinearGMRESSolver` is a wrapper class for the the Pysparse_
    `itsolvers.gmres()` and `precon.jacobi()` methods.

    .. _Pysparse: http://pysparse.sourceforge.net

    """

    solveFnc = staticmethod(krylov.gmres)

    DEFAULT_PRECONDITIONER = JacobiPreconditioner
