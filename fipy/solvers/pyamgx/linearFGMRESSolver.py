from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.preconditioners import JacobiPreconditioner

__all__ = ["LinearFGMRESSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearFGMRESSolver(PyAMGXSolver):
    """Interface to the Flexible Generalized Minimum RESidual
    (:term:`FGMRES`) solver in :ref:`PYAMGX`.

    Uses a
    :class:`~fipy.solvers.pyamgx.preconditioners.JacobiPreconditioner` by
    default.
    """

    CONFIG_DICT = {
        "config_version": 2,
        "determinism_flag": 1,
        "exception_handling" : 1,
        "solver": {
            "monitor_residual": 1,
            "solver": "FGMRES",
        }
    }

    DEFAULT_PRECONDITIONER = JacobiPreconditioner
