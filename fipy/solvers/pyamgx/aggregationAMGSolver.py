from __future__ import unicode_literals
from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.smoothers import BlockJacobiSmoother

__all__ = ["AggregationAMGSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class AggregationAMGSolver(PyAMGXSolver):
    """
    The `AggregationAMGSolver` is an interface to the aggregation AMG solver in
    AMGX, with a Jacobi smoother by default.
    """
    def __init__(self, tolerance=1e-10, iterations=2000,
                 precon=None,
                 smoother=BlockJacobiSmoother(),
                 **kwargs):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner to use.
          - `smoother`: Smoother to use.
          - `kwargs`: Keyword arguments specifying other AMGX solver options.
        """
        config_dict = {
            "config_version": 2,
            "determinism_flag": 1,
            "solver": {
                "algorithm": "AGGREGATION",
                "solver": "AMG",
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
                precon=precon,
                smoother=smoother,
                **kwargs)
