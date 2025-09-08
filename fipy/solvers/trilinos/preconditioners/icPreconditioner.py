__docformat__ = 'restructuredtext'

from PyTrilinos import IFPACK

from .trilinosPreconditioner import TrilinosPreconditioner

__all__ = ["ICPreconditioner"]

class ICPreconditioner(TrilinosPreconditioner):
    """Incomplete Cholesky Preconditioner from IFPACK for :class:`~fipy.solvers.trilinos.trilinosSolver.TrilinosSolver`.
    """

    def _applyToSolver(self, solver, matrix):
        Factory = IFPACK.Factory()
        Prec = Factory.Create("IC", matrix)
        Prec.Initialize()
        Prec.Compute()
        solver.SetPrecOperator(Prec)
