from fipy.solvers.pyamgx import PyAMGXSolver

__all__ = ["LinearPCGSolver"]

class LinearPCGSolver(PyAMGXSolver):

    def __init__(self, tolerance=1e-10, iterations=2000):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner to use.

        """
        config_dict = {
            "config_version": 2, 
            "determinism_flag": 1, 
            "solver": {
                "preconditioner": {
                    "print_grid_stats": 0, 
                    "algorithm": "AGGREGATION", 
                    "print_vis_data": 0, 
                    "solver": "AMG", 
                    "smoother": {
                        "relaxation_factor": 0.8, 
                        "scope": "jacobi", 
                        "solver": "BLOCK_JACOBI", 
                        "monitor_residual": 0, 
                        "print_solve_stats": 0
                    }, 
                    "print_solve_stats": 0, 
                    "presweeps": 0, 
                    "selector": "SIZE_2", 
                    "coarse_solver": "NOSOLVER", 
                    "max_iters": 1, 
                    "monitor_residual": 0, 
                    "store_res_history": 0, 
                    "scope": "amg", 
                    "max_levels": 100, 
                    "postsweeps": 3, 
                    "cycle": "V"
                }, 
                "solver": "PCG", 
                "print_solve_stats": 0, 
                "obtain_timings": 0, 
                "max_iters": 1000, 
                "monitor_residual": 1, 
                "convergence": "RELATIVE_INI", 
                "scope": "main", 
                "tolerance": 0.0001, 
                "norm": "L2"
            }
        }
        config_dict['solver']['tolerance'] = tolerance
        config_dict['solver']['max_iters'] = iterations

        super(LinearPCGSolver, self).__init__(config_dict, tolerance=tolerance, iterations=iterations)

    def _canSolveAsymmetric(self):
        return False
