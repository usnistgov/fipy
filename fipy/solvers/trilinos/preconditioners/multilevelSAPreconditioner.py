__docformat__ = 'restructuredtext'

from .multilevelPreconditioner import MultilevelPreconditioner

__all__ = ["MultilevelSAPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class MultilevelSAPreconditioner(MultilevelPreconditioner):
    """Classical smoothed aggregation multilevel preconditioner for :class:`~fipy.solvers.trilinos.trilinosSolver.TrilinosSolver`.

    Suitable for symmetric positive definite or nearly symmetric positive
    definite systems.
    """

    @property
    def _parameterList(self):
        return {
            "output": 0,
            "max levels": self.levels,
            "prec type": "MGV",
            "increasing or decreasing": "increasing",
            "aggregation: type": "Uncoupled-MIS",
            "aggregation: damping factor" : 4. / 3.,
            "eigen-analysis: type": "cg",
            "eigen-analysis: iterations": 10,
            "smoother: sweeps": 2,
            "smoother: damping factor": 1.0,
            "smoother: pre or post": "both",
            "smoother: type": "symmetric Gauss-Seidel",
            "coarse: type": "Amesos-KLU",
            "coarse: max size": 128
        }
