__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["DomDecompPreconditioner"]

class DomDecompPreconditioner(Preconditioner):
    """
    Domain Decomposition preconditioner for Trilinos solvers.

    """

    def _ApplyToSolver(self, solver, matrix):
        solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_dom_decomp)
