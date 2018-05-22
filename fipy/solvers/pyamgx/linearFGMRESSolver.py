from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.preconditioners import BlockJacobiPreconditioner

__all__ = ["LinearFGMRESSolver"]

class LinearFGMRESSolver(PyAMGXSolver):
    """
    The `LinearFGMRESSolver` is an interface to the FGMRES solver in
    AMGX, with a Jacobi preconditioner by default.
    """

    def __init__(self, tolerance=1e-10, iterations=2000,
                 preconditioner=BlockJacobiPreconditioner(),
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
                "monitor_residual": 1,
                "solver": "FGMRES",
                "preconditioner": {
                    "solver": "NOSOLVER"
                }
            }
        }
        super(LinearFGMRESSolver, self).__init__(
                config_dict,
                tolerance=tolerance,
                iterations=iterations,
                preconditioner=preconditioner,
                **kwargs)
