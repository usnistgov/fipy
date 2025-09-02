from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from .multilevelPreconditioner import MultilevelPreconditioner

__all__ = ["MultilevelDDMLPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class MultilevelDDMLPreconditioner(MultilevelPreconditioner):
    """3-level algebraic domain decomposition multilevel preconditioner for :class:`~fipy.solvers.trilinos.trilinosSolver.TrilinosSolver`.
    """

    def __init__(self, levels=3):
        """

        Parameters
        ----------
        levels : int
            Maximum number of levels
        """
        self.levels = levels

    @property
    def _parameterList(self):
        return {
            "output": 0,
            "max levels": self.levels,
            "prec type": "MGV",
            "increasing or decreasing": "increasing",
            "aggregation: type": "METIS",
            "aggregation: nodes per aggregate": 512,
            "aggregation: next-level aggregates per process": 128,
            "aggregation: damping factor" : 4. / 3.,
            "eigen-analysis: type": "power-method",
            "eigen-analysis: iterations": 20,
            "smoother: sweeps": 1,
            "smoother: pre or post": "both",
            "smoother: type": "Aztec",
            "smoother: Aztec as solver": False,
            "coarse: type": "Amesos-KLU",
            "coarse: max size": 128
        }
