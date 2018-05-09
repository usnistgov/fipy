from fipy.solvers.pyamgx import PyAMGXSolver

__all__ = ["LinearPCGSolver"]

class LinearPCGSolver(PyAMGXSolver):
    """
    The `LinearPCGSolver` is an interface to the PCG solver in
    AMGX, with no preconditioning by default.
    """

    def __init__(self, tolerance=1e-10, iterations=2000, preconditioner=None,
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
            "solver": {"solver": "PCG",
                       "preconditioner": "NOSOLVER"},
        }
        super(LinearPCGSolver, self).__init__(config_dict, tolerance=tolerance, iterations=iterations,
                preconditioner=preconditioner, **kwargs)

    def _canSolveAsymmetric(self):
        return False
