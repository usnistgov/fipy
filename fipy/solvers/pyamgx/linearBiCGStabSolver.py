from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.preconditioners import JacobiPreconditioner

__all__ = ["LinearBiCGStabSolver"]

class LinearBiCGStabSolver(PyAMGXSolver):
    """Interface to the Biconjugate Gradient (Stabilized) (:term:`BiCGSTAB`)
    solver in :ref:`PYAMGX`.

    Uses a
    :class:`~fipy.solvers.pyamgx.preconditioners.JacobiPreconditioner` by
    default.
    """

    CONFIG_DICT = {
        "config_version": 2,
        "determinism_flag": 1,
        "exception_handling" : 1,
        "solver": {
            "convergence": "RELATIVE_INI_CORE",
            "monitor_residual": 1,
            "solver": "PBICGSTAB",
        }
    }

    DEFAULT_PRECONDITIONER = JacobiPreconditioner
