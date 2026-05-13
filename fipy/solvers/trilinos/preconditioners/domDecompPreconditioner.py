__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from .trilinosPreconditioner import TrilinosPreconditioner

__all__ = ["DomDecompPreconditioner"]

class DomDecompPreconditioner(TrilinosPreconditioner):
    """Domain Decomposition preconditioner for :class:`~fipy.solvers.trilinos.trilinosSolver.TrilinosSolver`.
    """

    def _applyToSolver(self, solver, matrix):
        solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_dom_decomp)
