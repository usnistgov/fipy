


__docformat__ = 'restructuredtext'

from PyTrilinos import ML

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["MultilevelSolverSmootherPreconditioner"]

class MultilevelSolverSmootherPreconditioner(Preconditioner):
    """
    Multilevel preconditioner for Trilinos solvers using Aztec solvers
    as smoothers.

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
        self.Prec.SetParameterList({"output": 0, "smoother: type" : "Aztec", "smoother: Aztec as solver" : True})
        self.Prec.ComputePreconditioner()
        solver.SetPrecOperator(self.Prec)
