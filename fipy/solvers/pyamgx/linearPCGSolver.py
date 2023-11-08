from __future__ import unicode_literals
from fipy.solvers.pyamgx import PyAMGXSolver
from fipy.solvers.pyamgx.preconditioners import JacobiPreconditioner

__all__ = ["LinearPCGSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearPCGSolver(PyAMGXSolver):
    """
    The `LinearPCGSolver` is an interface to the PCG solver in
    AMGX, with no preconditioning by default.
    """

    def __init__(self, tolerance=1e-10, criterion="legacy",
                 iterations=1000, precon=JacobiPreconditioner(),
                 **kwargs):
        """
        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'unscaled', 'RHS', 'matrix', 'initial', 'legacy'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.pyamgx.preconditioners.Preconditioner, optional
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
                "solver": "PCG",
                "preconditioner": {
                   "solver": "NOSOLVER",
                }
            }
        }
        super(LinearPCGSolver, self).__init__(
                config_dict,
                tolerance=tolerance,
                criterion=criterion,
                iterations=iterations,
                precon=precon,
                **kwargs)

    def _canSolveAsymmetric(self):
        return False
