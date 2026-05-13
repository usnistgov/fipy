__docformat__ = 'restructuredtext'

from PyTrilinos import ML

from .multilevelPreconditioner import MultilevelPreconditioner

__all__ = ["MultilevelSolverSmootherPreconditioner"]

class MultilevelSolverSmootherPreconditioner(MultilevelPreconditioner):
    """Multilevel preconditioner using Aztec solvers as smoothers for :class:`~fipy.solvers.trilinos.trilinosSolver.TrilinosSolver`.
    """

    @property
    def _parameterList(self):
        return {
            "output": 0,
            "max levels": self.levels,
            "smoother: type": "Aztec",
            "smoother: Aztec as solver": True
        }
