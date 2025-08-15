from __future__ import unicode_literals
from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.preconditioners import JacobiPreconditioner

__all__ = ["LinearPCGSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearPCGSolver(PyAMGXSolver):
    """Interface to the preconditioned conjugate gradient (:term:`PCG`)
    solver in :ref:`PYAMGX`.

    Uses :class:`~fipy.solvers.pyamgx.preconditioners.JacobiPreconditioner` by default.
    """

    CONFIG_DICT = {
        "config_version": 2,
        "determinism_flag": 1,
        "exception_handling" : 1,
        "solver": {
            "convergence": "RELATIVE_INI_CORE",
            "monitor_residual": 1,
            "solver": "PCG",
            "preconditioner": {
               "solver": "NOSOLVER",
            }
        }
    }

    DEFAULT_PRECONDITIONER = JacobiPreconditioner

    def _canSolveAsymmetric(self):
        return False
