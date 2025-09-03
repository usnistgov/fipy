__docformat__ = 'restructuredtext'

from .multilevelPreconditioner import MultilevelPreconditioner

__all__ = ["MultilevelNSSAPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class MultilevelNSSAPreconditioner(MultilevelPreconditioner):
    """Energy-based minimizing smoothed aggregation preconditioner for :class:`~fipy.solvers.trilinos.trilinosSolver.TrilinosSolver`.

    Suitable for highly convective non-symmetric fluid flow problems.
    """

    @property
    def _parameterList(self):
        return {
            "output": 0,
            "max levels": self.levels,
            "prec type": "MGW",
            "increasing or decreasing": "increasing",
            "aggregation: type": "Uncoupled-MIS",
            "energy minimization: enable": True,
            "eigen-analysis: type": "power-method",
            "eigen-analysis: iterations": 20,
            "smoother: sweeps": 4,
            "smoother: damping factor": 0.67,
            "smoother: pre or post": "post",
            "smoother: type": "symmetric Gauss-Seidel",
            "coarse: type": "Amesos-KLU",
            "coarse: max size": 256
        }
