__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from .trilinosPreconditioner import TrilinosPreconditioner

__all__ = ["JacobiPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class JacobiPreconditioner(TrilinosPreconditioner):
    """Jacobi preconditioner for :class:`~fipy.solvers.trilinos.trilinosSolver.TrilinosSolver`.
    """

    def _applyToSolver(self, solver, matrix):
        solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_Jacobi)
