from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.preconditioners import JacobiPreconditioner

__all__ = ["LinearGMRESSolver"]

class LinearGMRESSolver(PyAMGXSolver):
    """Interface to the Generalized Minimum RESidual (:term:`GMRES`) solver in
    :ref:`PYAMGX`.

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
            "solver": "GMRES",
            "convergence": "RELATIVE_INI_CORE",
        }
    }

    DEFAULT_PRECONDITIONER = JacobiPreconditioner
