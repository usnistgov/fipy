from fipy.solvers.scipy.linearGMRESSolver import LinearGMRESSolver as ScipyLinearGMRESSolver
from fipy.solvers.pyAMG.preconditioners.smoothedAggregationPreconditioner import SmoothedAggregationPreconditioner

__all__ = ["LinearGMRESSolver"]

class LinearGMRESSolver(ScipyLinearGMRESSolver):
    """
    The `LinearGMRESSolver` is an interface to the GMRES solver in
    Scipy, using the pyAMG `SmoothedAggregationPreconditioner` by
    default.
    """

    def __init__(self, tolerance=1e-15, iterations=2000, precon=SmoothedAggregationPreconditioner()):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner to use.

        """

        super(LinearGMRESSolver, self).__init__(tolerance=tolerance, iterations=iterations, precon=precon)
