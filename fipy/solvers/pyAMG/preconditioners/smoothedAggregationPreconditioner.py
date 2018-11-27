


from pyamg import smoothed_aggregation_solver

__all__ = ["SmoothedAggregationPreconditioner"]

class SmoothedAggregationPreconditioner():
    def __init__(self):
        pass

    def _applyToMatrix(self, A):
        return smoothed_aggregation_solver(A).aspreconditioner(cycle='V')
