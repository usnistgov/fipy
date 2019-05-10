from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["JacobiPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class JacobiPreconditioner(Preconditioner):
    """
    Jacobi Preconditioner for Trilinos solvers.

    """

    def _applyToSolver(self, solver, matrix):
        solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_Jacobi)
