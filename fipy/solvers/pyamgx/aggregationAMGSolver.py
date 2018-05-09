from fipy.solvers.pyamgx import PyAMGXSolver

__all__ = ["AggregationAMGSolver"]

class AggregationAMGSolver(PyAMGXSolver):
    """
    The `AggregationAMGSolver` is an interface to the classical AMG solver in
    AMGX, with a Gauss Siedel smoother by default.
    """
    def __init__(self, tolerance=1e-10, iterations=1000,
                 preconditioner=None,
                 smoother=None,
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
                "smoother": "MULTICOLOR_GS", 
                "coarse_solver": "NOSOLVER",
                "symmetric_GS": 1, 
                "selector": "SIZE_2", 
                "monitor_residual": 1, 
                "max_levels": 1000, 
                "matrix_coloring_scheme": "MIN_MAX", 
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
