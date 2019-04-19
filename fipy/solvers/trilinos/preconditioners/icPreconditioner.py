__docformat__ = 'restructuredtext'

from PyTrilinos import IFPACK

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["ICPreconditioner"]

class ICPreconditioner(Preconditioner):
    """
    Incomplete Cholesky Preconditioner from IFPACK for Trilinos Solvers.

    """

    def _applyToSolver(self, solver, matrix):
        Factory = IFPACK.Factory()
        Prec = Factory.Create("IC", matrix)
        Prec.Initialize()
        Prec.Compute()
        solver.SetPrecOperator(Prec)
