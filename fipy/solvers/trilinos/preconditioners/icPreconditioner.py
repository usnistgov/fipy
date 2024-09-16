from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import IFPACK

from .trilinosPreconditioner import TrilinosPreconditioner

__all__ = ["ICPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class ICPreconditioner(TrilinosPreconditioner):
    """Incomplete Cholesky Preconditioner from IFPACK for :class:`~fipy.solvers.trilinos.trilinosSolver.TrilinosSolver`.
    """

    def _applyToSolver(self, solver, matrix):
        Factory = IFPACK.Factory()
        Prec = Factory.Create(text_to_native_str("IC"), matrix)
        Prec.Initialize()
        Prec.Compute()
        solver.SetPrecOperator(Prec)
