__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from fipy.solvers.trilinos.trilinosAztecOOSolver import TrilinosAztecOOSolver
from fipy.solvers.trilinos.preconditioners.multilevelDDPreconditioner import MultilevelDDPreconditioner

__all__ = ["LinearGMRESSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearGMRESSolver(TrilinosAztecOOSolver):

    """Interface to the generalized minimal residual (:term:`GMRES`) solver
    in :ref:`TRILINOS`.

    Uses the
    :class:`~fipy.solvers.trilinos.preconditioners.multilevelDDPreconditioner.MultilevelDDPreconditioner`
    by default.
    """

    solver = AztecOO.AZ_gmres

    DEFAULT_PRECONDITIONER = MultilevelDDPreconditioner
