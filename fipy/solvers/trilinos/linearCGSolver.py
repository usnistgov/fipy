from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from fipy.solvers.trilinos.trilinosAztecOOSolver import TrilinosAztecOOSolver
from fipy.solvers.trilinos.preconditioners.multilevelDDPreconditioner import MultilevelDDPreconditioner

__all__ = ["LinearCGSolver", "LinearPCGSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearCGSolver(TrilinosAztecOOSolver):

    """Interface to the conjugate gradient (:term:`CG`) solver in
    :ref:`TRILINOS`.

    Uses the
    :class:`~fipy.solvers.trilinos.preconditioners.multilevelDDPreconditioner.MultilevelDDPreconditioner`
    by default.
    """

    solver = AztecOO.AZ_cg

    DEFAULT_PRECONDITIONER = MultilevelDDPreconditioner

    def _canSolveAsymmetric(self):
        return False

LinearPCGSolver = LinearCGSolver
