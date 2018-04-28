from fipy.solvers.pyamgx import PyAMGXSolver

__all__ = ["LinearGMRESSolver"]

class LinearGMRESSolver(PyAMGXSolver):

    def __init__(self, tolerance=1e-10, iterations=2000):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner to use.

        """
        cfgDict = {
            "config_version": 2, 
            "determinism_flag": 1,
            "exception_handling" : 1,
            "solver": {
                "store_res_history": 1, 
                "solver": "GMRES", 
                "obtain_timings": 1, 
                "preconditioner": {
                    "interpolator": "D2", 
                    "solver": "AMG", 
                    "smoother": "JACOBI_L1", 
                    "presweeps": 2, 
                    "selector": "PMIS", 
                    "coarsest_sweeps": 2, 
                    "coarse_solver": "NOSOLVER", 
                    "max_iters": 1, 
                    "interp_max_elements": 4, 
                    "min_coarse_rows": 2, 
                    "scope": "amg_solver", 
                    "max_levels": 24, 
                    "cycle": "V", 
                    "postsweeps": 2
                }, 
                "max_iters": 100, 
                "monitor_residual": 1, 
                "gmres_n_restart": 10, 
                "convergence": "RELATIVE_INI_CORE", 
                "tolerance": 1e-06, 
                "norm": "L2"
            }
        }
        cfgDict['solver']['tolerance'] = tolerance
        cfgDict['solver']['max_iters'] = iterations

        super(LinearGMRESSolver, self).__init__(cfgDict)
