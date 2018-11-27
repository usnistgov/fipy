


__docformat__ = 'restructuredtext'

from PyTrilinos import ML

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["MultilevelNSSAPreconditioner"]

class MultilevelNSSAPreconditioner(Preconditioner):
    """
    Energy-based minimizing smoothed aggregation suitable for highly
    convective non-symmetric fluid flow problems.
    """
    def _applyToSolver(self, solver, matrix):
        if matrix.NumGlobalNonzeros() <= matrix.NumGlobalRows():
            return

        self.Prec = ML.MultiLevelPreconditioner(matrix, False)

        self.Prec.SetParameterList({"output": 0,
                                    "max levels" : 10,
                                    "prec type" : "MGW",
                                    "increasing or decreasing" : "increasing",
                                    "aggregation: type" : "Uncoupled-MIS",
                                    "energy minimization: enable" : True,
                                    "eigen-analysis: type" : "power-method",
                                    "eigen-analysis: iterations" : 20,
                                    "smoother: sweeps" : 4,
                                    "smoother: damping factor" : 0.67,
                                    "smoother: pre or post" : 'post',
                                    "smoother: type" : "symmetric Gauss-Seidel",
                                    "coarse: type" : 'Amesos-KLU',
                                    "coarse: max size" : 256
                                    })

        self.Prec.ComputePreconditioner()

        solver.SetPrecOperator(self.Prec)
