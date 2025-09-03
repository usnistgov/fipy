__docformat__ = 'restructuredtext'

from .multilevelPreconditioner import MultilevelPreconditioner

__all__ = ["MultilevelSGSPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class MultilevelSGSPreconditioner(MultilevelPreconditioner):
    """Multilevel preconditioner using Symmetric Gauss-Seidel smoothing for :class:`~fipy.solvers.trilinos.trilinosSolver.TrilinosSolver`.
    """

    @property
    def _parameterList(self):
        return {
            "output": 0,
            "max levels": self.levels,
            "smoother: type": "symmetric Gauss-Seidel",
        }
