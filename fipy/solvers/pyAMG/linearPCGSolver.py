#!/usr/bin/env python


from fipy.solvers.scipy.linearPCGSolver import LinearPCGSolver as ScipyLinearPCGSolver
from fipy.solvers.pyAMG.preconditioners.smoothedAggregationPreconditioner import SmoothedAggregationPreconditioner

__all__ = ["LinearPCGSolver"]

class LinearPCGSolver(ScipyLinearPCGSolver):
    """
    The `LinearPCGSolver` is an interface to the PCG solver in Scipy,
    using the pyAMG `SmoothedAggregationPreconditioner` by default.
    """

    def __init__(self, tolerance=1e-15, iterations=2000, precon=SmoothedAggregationPreconditioner()):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner to use.

        """

        super(LinearPCGSolver, self).__init__(tolerance=tolerance, iterations=iterations, precon=precon)
