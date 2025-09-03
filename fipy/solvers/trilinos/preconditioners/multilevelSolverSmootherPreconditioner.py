__docformat__ = 'restructuredtext'

from PyTrilinos import ML

from .multilevelPreconditioner import MultilevelPreconditioner

__all__ = ["MultilevelSolverSmootherPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

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
