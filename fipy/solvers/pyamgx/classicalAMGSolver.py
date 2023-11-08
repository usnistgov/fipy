from __future__ import unicode_literals
from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.smoothers import BlockJacobiSmoother

__all__ = ["ClassicalAMGSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class ClassicalAMGSolver(PyAMGXSolver):
    """
    The `ClassicalAMGSolver` is an interface to the classical AMG solver in
    AMGX, with a Jacobi smoother by default.
    """
    def __init__(self, tolerance=1e-10, criterion="default",
                 iterations=1000, precon=None,
                 smoother=BlockJacobiSmoother(),
                 **kwargs):
        """
        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'unscaled', 'RHS', 'matrix', 'initial', 'legacy'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.pyamgx.preconditioners.Preconditioner, optional
        smoother : ~fipy.solvers.pyamgx.smoothers.Smoother, optional
        **kwargs
            Other AMGX solver options
        """
        config_dict = {
            "config_version": 2,
            "determinism_flag": 1,
            "solver": {
                "algorithm": "CLASSICAL",
                "solver": "AMG",
                "monitor_residual": 1,
                "max_levels": 1000,
            }
        }
        super(ClassicalAMGSolver, self).__init__(
                config_dict,
                tolerance=tolerance,
                criterion=criterion,
                iterations=iterations,
                precon=precon,
                smoother=smoother,
                **kwargs)
