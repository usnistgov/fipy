from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import ML

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["MultilevelSGSPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class MultilevelSGSPreconditioner(Preconditioner):
    """
    Multilevel preconditioner for Trilinos solvers using Symmetric Gauss-Seidel smoothing.

    """
    def __init__(self, levels=10):
        """
        Initialize the multilevel preconditioner

        - `levels`: Maximum number of levels
        """
        self.levels = levels

    def _applyToSolver(self, solver, matrix):
        if matrix.NumGlobalNonzeros() <= matrix.NumGlobalRows():
            return

        self.Prec = ML.MultiLevelPreconditioner(matrix, False)
        self.Prec.SetParameterList({text_to_native_str("output"): 0, text_to_native_str("smoother: type") : text_to_native_str("symmetric Gauss-Seidel")})
        self.Prec.ComputePreconditioner()
        solver.SetPrecOperator(self.Prec)
