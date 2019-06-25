from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["DomDecompPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class DomDecompPreconditioner(Preconditioner):
    """
    Domain Decomposition preconditioner for Trilinos solvers.

    """

    def _ApplyToSolver(self, solver, matrix):
        solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_dom_decomp)
