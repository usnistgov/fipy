from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from .trilinosPreconditioner import TrilinosPreconditioner

__all__ = ["DomDecompPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class DomDecompPreconditioner(TrilinosPreconditioner):
    """Domain Decomposition preconditioner for :class:`~fipy.solvers.trilinos.trilinosSolver.TrilinosSolver`.
    """

    def _applyToSolver(self, solver, matrix):
        solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_dom_decomp)
