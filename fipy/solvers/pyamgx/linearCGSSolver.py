from fipy.solvers.pyamgx import PyAMGXSolver

__all__ = ["LinearCGSSolver"]

class LinearCGSSolver(PyAMGXSolver):

    def __init__(self, tolerance=1e-10, iterations=2000):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner to use.

        """
        cfgDict = {
            "config_version": 2, 
            "solver": {
                "preconditioner": {
                    "error_scaling": 0, 
                    "algorithm": "AGGREGATION", 
                    "solver": "AMG", 
                    "smoother": "BLOCK_JACOBI", 
                    "presweeps": 0, 
                    "selector": "SIZE_2", 
                    "coarse_solver": "NOSOLVER", 
                    "max_iters": 1, 
                    "min_coarse_rows": 32, 
                    "relaxation_factor": 0.75, 
                    "scope": "amg", 
                    "max_levels": 100, 
                    "postsweeps": 3, 
                    "cycle": "V"
                }, 
                "use_scalar_norm": 1, 
                "solver": "PCG", 
                "obtain_timings": 1, 
                "max_iters": 100, 
                "monitor_residual": 1, 
                "gmres_n_restart": 10, 
                "convergence": "RELATIVE_INI_CORE", 
                "scope": "main", 
                "tolerance": 1e-10, 
                "norm": "L2"
            }
        }
        cfgDict['solver']['tolerance'] = tolerance
        cfgDict['solver']['max_iters'] = iterations

        super(LinearCGSSolver, self).__init__(cfgDict)
