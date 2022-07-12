from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from .domDecompPreconditioner import DomDecompPreconditioner

__all__ = ["ILUPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class ILUPreconditioner(DomDecompPreconditioner):
    """
    ILU Domain Decomposition preconditioner for Trilinos solvers.

    """

    def _applyToSolver(self, solver, matrix):
        super(ILUPreconditioner, self)._applyToSolver(solver, matrix)
        solver.SetAztecOption(AztecOO.AZ_subdomain_solve, AztecOO.AZ_ilu)
