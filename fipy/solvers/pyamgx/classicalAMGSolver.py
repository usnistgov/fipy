from fipy.solvers.pyamgx import PyAMGXSolver

__all__ = ["ClassicalAMGSolver"]

class ClassicalAMGSolver(PyAMGXSolver):
    """
    The `ClassicalAMGSolver` is an interface to the classical AMG solver in
    AMGX, with a Jacobi smoother by default.
    """
    def __init__(self, tolerance=1e-10, iterations=2000,
                 preconditioner=None,
                 smoother=None,
                 **kwargs):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `preconditioner`: Preconditioner to use.
          - `smoother`: Smoother to use.
          - `kwargs`: Keyword arguments specifying other AMGX solver options.
        """
        config_dict = {
            "config_version": 2, 
            "determinism_flag": 1, 
            "solver": {
                "algorithm": "CLASSICAL", 
                "solver": "AMG", 
                "coarse_solver": "NOSOLVER",
                "monitor_residual": 1, 
                "print_solve_stats": 1, 
                "max_levels": 1000, 
            }
        }
        super(ClassicalAMGSolver, self).__init__(
                config_dict,
                tolerance=tolerance,
                iterations=iterations,
                preconditioner=preconditioner,
                smoother=smoother,
                **kwargs)
