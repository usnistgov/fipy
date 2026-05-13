from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.smoothers import BlockJacobiSmoother

__all__ = ["ClassicalAMGSolver"]

class ClassicalAMGSolver(PyAMGXSolver):
    """Interface to the classical algebraic multigrid (:term:`AMG`) solver
    in :ref:`PYAMGX`.

    Uses a :class:`~fipy.solvers.pyamgx.smoothers.BlockJacobiSmoother`
    smoother by default.
    """

    CONFIG_DICT = {
        "config_version": 2,
        "determinism_flag": 1,
        "solver": {
            "algorithm": "CLASSICAL",
            "solver": "AMG",
            "monitor_residual": 1,
            "max_levels": 1000,
        }
    }

    DEFAULT_SMOOTHER = BlockJacobiSmoother
