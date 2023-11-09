from __future__ import unicode_literals
from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.preconditioners import JacobiPreconditioner

__all__ = ["LinearGMRESSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearGMRESSolver(PyAMGXSolver):
    """
    The `LinearGMRESSolver` is an interface to the GMRES solver in
    AMGX, with a Jacobi preconditioner by default.
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
