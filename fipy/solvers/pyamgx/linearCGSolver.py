from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.preconditioners import *

__all__ = ["LinearCGSolver", "LinearPCGSolver"]

class LinearCGSolver(PyAMGXSolver):
    """
    The `LinearCGSolver` is an interface to the PCG solver in
    AMGX, with no preconditioning by default.
    """

    def __init__(self, tolerance=1e-10, iterations=2000, preconditioner=BlockJacobiPreconditioner(),
            **kwargs):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `preconditioner`: Preconditioner to use.
          - `kwargs`: Keyword arguments specifying other AMGX solver options.
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
                preconditioner=preconditioner,
                **kwargs)

    def _canSolveAsymmetric(self):
        return False

LinearPCGSolver = LinearCGSolver
