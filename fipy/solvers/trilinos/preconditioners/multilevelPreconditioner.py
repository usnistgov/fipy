__docformat__ = 'restructuredtext'

from PyTrilinos import ML

from .trilinosPreconditioner import TrilinosPreconditioner

__all__ = ["MultilevelPreconditioner"]

class MultilevelPreconditioner(TrilinosPreconditioner):
    """Base class for multilevel preconditioners  for :class:`~fipy.solvers.trilinos.trilinosSolver.TrilinosSolver`.
    """

    def __init__(self, levels=10):
        """
        Parameters
        ----------
        levels : int
            Maximum number of levels
        """
        self.levels = levels

    @property
    def _parameterList(self):
        """Trilinos preconditioner parameters.

        Implemented as a property to avoid
        `side-effects <https://docs.python-guide.org/writing/gotchas/#mutable-default-arguments>`_.

        Returns
        -------
        dict
        """
        raise NotImplementedError

    def _applyToSolver(self, solver, matrix):
        if matrix.NumGlobalNonzeros() <= matrix.NumGlobalRows():
            return

        self.Prec = ML.MultiLevelPreconditioner(matrix, False)

        self.Prec.SetParameterList(self._parameterList)

        self.Prec.ComputePreconditioner()

        solver.SetPrecOperator(self.Prec)
