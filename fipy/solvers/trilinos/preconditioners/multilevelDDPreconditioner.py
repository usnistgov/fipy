__docformat__ = 'restructuredtext'

from .multilevelPreconditioner import MultilevelPreconditioner

__all__ = ["MultilevelDDPreconditioner"]

class MultilevelDDPreconditioner(MultilevelPreconditioner):
    """Classical smoothed aggregation-based 2-level domain decomposition preconditioner for :class:`~fipy.solvers.trilinos.trilinosSolver.TrilinosSolver`.
    """

    def __init__(self, levels=2):
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
            "aggregation: local aggregates": 1,
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
