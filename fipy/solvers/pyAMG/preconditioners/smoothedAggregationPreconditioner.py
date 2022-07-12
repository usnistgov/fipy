from __future__ import unicode_literals
from builtins import object

from pyamg import smoothed_aggregation_solver

from ...scipy.preconditioners.preconditioner import Preconditioner

__all__ = ["SmoothedAggregationPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class SmoothedAggregationPreconditioner(Preconditioner):

    def _applyToMatrix(self, A):
        return smoothed_aggregation_solver(A).aspreconditioner(cycle='V')
