from __future__ import unicode_literals
from builtins import object
from pyamg import smoothed_aggregation_solver

__all__ = ["SmoothedAggregationPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class SmoothedAggregationPreconditioner(object):
    """Preconditioner based on `PyAMG smoothed_aggregation_solver <https://pyamg.readthedocs.io/en/latest/generated/pyamg.aggregation.html#pyamg.aggregation.smoothed_aggregation_solver>`_
    """
    def __init__(self):
        pass

    def _applyToMatrix(self, A):
        return smoothed_aggregation_solver(A).aspreconditioner(cycle='V')
