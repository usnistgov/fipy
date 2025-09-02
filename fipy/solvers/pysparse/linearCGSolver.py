from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from pysparse.itsolvers import krylov

from .linearRHSSolver import LinearRHSSolver
from .preconditioners import SSORPreconditioner

__all__ = ["LinearCGSolver", "LinearPCGSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearCGSolver(LinearRHSSolver):
    """Interface to conjugate gradient method (:term:`CG`)
    of :ref:`Pysparse`.

    Uses :class:`~fipy.solvers.pysparse.preconditioners.SSORPreconditioner`
    by default.
    """

    solveFnc = staticmethod(krylov.pcg)

    DEFAULT_PRECONDITIONER = SSORPreconditioner

    def _canSolveAsymmetric(self):
        return False

LinearPCGSolver = LinearCGSolver
