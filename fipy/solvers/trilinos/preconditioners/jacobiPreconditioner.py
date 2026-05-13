__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from .trilinosPreconditioner import TrilinosPreconditioner

__all__ = ["JacobiPreconditioner"]

class JacobiPreconditioner(TrilinosPreconditioner):
    """Jacobi preconditioner for :class:`~fipy.solvers.trilinos.trilinosSolver.TrilinosSolver`.
    """

    def _applyToSolver(self, solver, matrix):
        solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_Jacobi)
