from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from pysparse.itsolvers import krylov

from .linearInitialSolver import LinearInitialSolver
from .preconditioners import JacobiPreconditioner

__all__ = ["LinearGMRESSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearGMRESSolver(LinearInitialSolver):
    """Interface to the Generalized Minimum RESidual (:term:`GMRES`) solver
    in :ref:`Pysparse`.

    Uses :class:`~fipy.solvers.pysparse.preconditioners.JacobiPreconditioner` by default.
    """

    solveFnc = staticmethod(krylov.gmres)

    DEFAULT_PRECONDITIONER = JacobiPreconditioner
