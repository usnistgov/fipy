from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import ML

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["MultilevelDDMLPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class MultilevelDDMLPreconditioner(Preconditioner):
    """
    Multilevel preconditioner for Trilinos solvers. 3-level algebraic domain decomposition.
    """

    def _applyToSolver(self, solver, matrix):
        if matrix.NumGlobalNonzeros() <= matrix.NumGlobalRows():
            return

        self.Prec = ML.MultiLevelPreconditioner(matrix, False)

        self.Prec.SetParameterList({text_to_native_str("output"): 0,
                                    text_to_native_str("max levels") : 3,
                                    text_to_native_str("prec type") : text_to_native_str("MGV"),
                                    text_to_native_str("increasing or decreasing") : text_to_native_str("increasing"),
                                    text_to_native_str("aggregation: type") : text_to_native_str("METIS"),
                                    text_to_native_str("aggregation: nodes per aggregate") : 512,
                                    text_to_native_str("aggregation: next-level aggregates per process") : 128,
                                    text_to_native_str("aggregation: damping factor") : 4. / 3.,
                                    text_to_native_str("eigen-analysis: type") : text_to_native_str("power-method"),
                                    text_to_native_str("eigen-analysis: iterations") : 20,
                                    text_to_native_str("smoother: sweeps") : 1,
                                    text_to_native_str("smoother: pre or post") : text_to_native_str("both"),
                                    text_to_native_str("smoother: type") : text_to_native_str("Aztec"),
                                    text_to_native_str("smoother: Aztec as solver") : False,
                                    text_to_native_str("coarse: type") : text_to_native_str("Amesos-KLU"),
                                    text_to_native_str("coarse: max size") : 128
                                    })

        self.Prec.ComputePreconditioner()

        solver.SetPrecOperator(self.Prec)
