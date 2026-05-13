__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from .domDecompPreconditioner import DomDecompPreconditioner

__all__ = ["ILUPreconditioner"]

class ILUPreconditioner(DomDecompPreconditioner):
    """Incomplete LU Domain Decomposition preconditioner for :class:`~fipy.solvers.trilinos.trilinosSolver.TrilinosSolver`.
    """

    def _applyToSolver(self, solver, matrix):
        super(ILUPreconditioner, self)._applyToSolver(solver, matrix)
        solver.SetAztecOption(AztecOO.AZ_subdomain_solve, AztecOO.AZ_ilu)
