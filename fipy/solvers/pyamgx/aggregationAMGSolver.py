from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.smoothers import BlockJacobiSmoother

__all__ = ["AggregationAMGSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class AggregationAMGSolver(PyAMGXSolver):
    """Interface to the aggregation algebraic multigrid (:term:`AMG`) solver
    in :ref:`PYAMGX`.

    Uses a :class:`~fipy.solvers.pyamgx.smoothers.BlockJacobiSmoother`
    smoother by default.
    """

    CONFIG_DICT = {
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

    DEFAULT_SMOOTHER = BlockJacobiSmoother
