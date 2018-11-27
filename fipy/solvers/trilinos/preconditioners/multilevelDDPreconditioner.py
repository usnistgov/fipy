


__docformat__ = 'restructuredtext'

from PyTrilinos import ML

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["MultilevelDDPreconditioner"]

class MultilevelDDPreconditioner(Preconditioner):
    """
    Multilevel preconditioner for Trilinos solvers. A classical smoothed
    aggregation-based 2-level domain decomposition.
    """

    def _applyToSolver(self, solver, matrix):
        if matrix.NumGlobalNonzeros() <= matrix.NumGlobalRows():
            return

        self.Prec = ML.MultiLevelPreconditioner(matrix, False)

        self.Prec.SetParameterList({"output": 0,
                                    "max levels" : 2,
                                    "prec type" : "MGV",
                                    "increasing or decreasing" : "increasing",
                                    "aggregation: type" : "METIS",
                                    "aggregation: local aggregates" : 1,
                                    "aggregation: damping factor" : 4. / 3.,
                                    "eigen-analysis: type" : "power-method",
                                    "eigen-analysis: iterations" : 20,
                                    "smoother: sweeps" : 1,
                                    "smoother: pre or post" : 'both',
                                    "smoother: type" : "Aztec",
                                    "smoother: Aztec as solver" : False,
                                    "coarse: type" : 'Amesos-KLU',
                                    "coarse: max size" : 128
                                    })

        self.Prec.ComputePreconditioner()

        solver.SetPrecOperator(self.Prec)
