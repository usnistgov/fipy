
from pyamg import smoothed_aggregation_solver

from ...scipy.preconditioners.scipyPreconditioner import ScipyPreconditioner

__all__ = ["SmoothedAggregationPreconditioner"]

class SmoothedAggregationPreconditioner(ScipyPreconditioner):
    """Preconditioner based on `PyAMG smoothed_aggregation_solver`_ for :class:`~fipy.solvers.scipy.scipySolver.ScipySolver`.

    .. _PyAMG smoothed_aggregation_solver: https://pyamg.readthedocs.io/en/latest/generated/pyamg.aggregation.html#pyamg.aggregation.smoothed_aggregation_solver
    """

    def __init__(self):
        pass

    def _applyToMatrix(self, matrix):
        return smoothed_aggregation_solver(matrix).aspreconditioner(cycle='V'), None
