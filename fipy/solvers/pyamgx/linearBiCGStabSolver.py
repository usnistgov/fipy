from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.preconditioners import *

__all__ = ["LinearBiCGStabSolver"]

class LinearBiCGStabSolver(PyAMGXSolver):
    """
    The `LinearBiCGStabSolver` is an interface to the PBICGSTAB solver in
    AMGX, with a Jacobi preconditioner by default.
    """

    def __init__(self, tolerance=1e-10, iterations=2000,
            precon=BlockJacobiPreconditioner(),
                 **kwargs):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner to use.
          - `kwargs`: Keyword arguments specifying other AMGX solver options.
        """
        config_dict = {
            "config_version": 2,
            "determinism_flag": 1,
            "exception_handling" : 1,
            "solver": {
                "convergence": "RELATIVE_INI_CORE",
                "monitor_residual": 1,
                "solver": "PBICGSTAB",
                "preconditioner": {
                    "solver": "NOSOLVER"
                }
            }
        }
        super(LinearBiCGStabSolver, self).__init__(
                config_dict,
                tolerance=tolerance,
                iterations=iterations,
                precon=precon,
                **kwargs)
