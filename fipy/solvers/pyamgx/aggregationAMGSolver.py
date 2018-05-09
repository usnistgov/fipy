from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.smoothers import BlockJacobiSmoother

__all__ = ["AggregationAMGSolver"]

class AggregationAMGSolver(PyAMGXSolver):
    """
    The `AggregationAMGSolver` is an interface to the aggregation AMG solver in
    AMGX, with a Jacobi smoother by default.
    """
    def __init__(self, tolerance=1e-10, iterations=2000,
                 preconditioner=None,
                 smoother=BlockJacobiSmoother(),
                 **kwargs):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `preconditioner`: Preconditioner to use.
          - `smoother`: Smoother to use.
          - `kwargs`: Keyword arguments specifying other AMGX solver options.
        """
        config_dict = {
            "config_version": 2, 
            "determinism_flag": 1, 
            "solver": {
                "algorithm": "AGGREGATION", 
                "solver": "AMG", 
                "coarse_solver": "NOSOLVER",
                "selector": "SIZE_2", 
                "monitor_residual": 1, 
                "max_levels": 1000, 
                "cycle": "V"
            }
        }
        super(AggregationAMGSolver, self).__init__(
                config_dict,
                tolerance=tolerance,
                iterations=iterations,
                preconditioner=preconditioner,
                smoother=smoother,
                **kwargs)
