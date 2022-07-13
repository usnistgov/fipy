from __future__ import unicode_literals
from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.preconditioners import JacobiPreconditioner

__all__ = ["LinearBiCGStabSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearBiCGStabSolver(PyAMGXSolver):
    """
    The `LinearBiCGStabSolver` is an interface to the PBICGSTAB solver in
    AMGX, with a Jacobi preconditioner by default.
    """

    def __init__(self, tolerance=1e-10, criterion="default",
                 iterations=1000, precon=JacobiPreconditioner(),
                 **kwargs):
        """
        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'unscaled', 'RHS', 'matrix', 'initial'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
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
                "solver": "PBICGSTAB",
                "preconditioner": {
                    "solver": "NOSOLVER"
                }
            }
        }
        super(LinearBiCGStabSolver, self).__init__(
                config_dict,
                criterion=criterion,
                tolerance=tolerance,
                iterations=iterations,
                precon=precon,
                **kwargs)
