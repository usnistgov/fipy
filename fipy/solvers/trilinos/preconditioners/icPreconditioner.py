from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import IFPACK

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["ICPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class ICPreconditioner(Preconditioner):
    """
    Incomplete Cholesky Preconditioner from IFPACK for Trilinos Solvers.

    """

    def _applyToSolver(self, solver, matrix):
        Factory = IFPACK.Factory()
        Prec = Factory.Create(text_to_native_str("IC"), matrix)
        Prec.Initialize()
        Prec.Compute()
        solver.SetPrecOperator(Prec)
