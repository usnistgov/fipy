


__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["JacobiPreconditioner"]

class JacobiPreconditioner(Preconditioner):
    """
    Jacobi Preconditioner for Trilinos solvers.

    """

    def _applyToSolver(self, solver, matrix):
        solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_Jacobi)
