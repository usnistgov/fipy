from __future__ import unicode_literals
from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.preconditioners import *

__all__ = ["LinearCGSolver", "LinearPCGSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearCGSolver(PyAMGXSolver):
    """
    The `LinearCGSolver` is an interface to the PCG solver in
    AMGX, with no preconditioning by default.
    """

    def __init__(self, tolerance=1e-10, iterations=2000, precon=BlockJacobiPreconditioner(),
            **kwargs):
        """
        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.pyamgx.preconditioners.preconditioners.Preconditioner, optional
        **kwargs
            Other AMGX solver options
        """
        config_dict = {
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
        super(LinearCGSolver, self).__init__(
                config_dict,
                tolerance=tolerance,
                iterations=iterations,
                precon=precon,
                **kwargs)

    def _canSolveAsymmetric(self):
        return False

LinearPCGSolver = LinearCGSolver
