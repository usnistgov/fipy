from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import ML

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["MultilevelNSSAPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class MultilevelNSSAPreconditioner(Preconditioner):
    """
    Energy-based minimizing smoothed aggregation suitable for highly
    convective non-symmetric fluid flow problems.
    """
    def _applyToSolver(self, solver, matrix):
        if matrix.NumGlobalNonzeros() <= matrix.NumGlobalRows():
            return

        self.Prec = ML.MultiLevelPreconditioner(matrix, False)

        self.Prec.SetParameterList({text_to_native_str("output"): 0,
                                    text_to_native_str("max levels") : 10,
                                    text_to_native_str("prec type") : text_to_native_str("MGW"),
                                    text_to_native_str("increasing or decreasing") : text_to_native_str("increasing"),
                                    text_to_native_str("aggregation: type") : text_to_native_str("Uncoupled-MIS"),
                                    text_to_native_str("energy minimization: enable") : True,
                                    text_to_native_str("eigen-analysis: type") : text_to_native_str("power-method"),
                                    text_to_native_str("eigen-analysis: iterations") : 20,
                                    text_to_native_str("smoother: sweeps") : 4,
                                    text_to_native_str("smoother: damping factor") : 0.67,
                                    text_to_native_str("smoother: pre or post") : text_to_native_str("post"),
                                    text_to_native_str("smoother: type") : text_to_native_str("symmetric Gauss-Seidel"),
                                    text_to_native_str("coarse: type") : text_to_native_str("Amesos-KLU"),
                                    text_to_native_str("coarse: max size") : 256
                                    })

        self.Prec.ComputePreconditioner()

        solver.SetPrecOperator(self.Prec)
