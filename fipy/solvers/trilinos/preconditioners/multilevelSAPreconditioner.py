#!/usr/bin/env python


__docformat__ = 'restructuredtext'

from PyTrilinos import ML

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["MultilevelSAPreconditioner"]

class MultilevelSAPreconditioner(Preconditioner):
    """
    Multilevel preconditioner for Trilinos solvers suitable classical
    smoothed aggregation for symmetric positive definite or nearly
    symmetric positive definite systems.
    """

    def _applyToSolver(self, solver, matrix):
        if matrix.NumGlobalNonzeros() <= matrix.NumGlobalRows():
            return

        self.Prec = ML.MultiLevelPreconditioner(matrix, False)

        self.Prec.SetParameterList({"output": 0,
                                    "max levels" : 10,
                                    "prec type" : "MGV",
                                    "increasing or decreasing" : "increasing",
                                    "aggregation: type" : "Uncoupled-MIS",
                                    "aggregation: damping factor" : 4. / 3.,
##                                    "energy minimization: enable" : False,
##                                    "smoother: type" : "Aztec",
##                                    "smoother: type" : "symmetric Gauss-Seidel",
##                                    "eigen-analysis: type" : "power-method",
                                    "eigen-analysis: type" : "cg",
                                    "eigen-analysis: iterations" : 10,
                                    "smoother: sweeps" : 2,
                                    "smoother: damping factor" : 1.0,
                                    "smoother: pre or post" : 'both',
                                    "smoother: type" : "symmetric Gauss-Seidel",
                                    "coarse: type" : 'Amesos-KLU',
                                    "coarse: max size" : 128
                                    })

        self.Prec.ComputePreconditioner()

        solver.SetPrecOperator(self.Prec)
