from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from .trilinosSolver import TrilinosSolver
from .preconditioners.jacobiPreconditioner import JacobiPreconditioner

__all__ = ["TrilinosAztecOOSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class TrilinosAztecOOSolver(TrilinosSolver):

    """
    .. attention:: This class is abstract, always create on of its subclasses.
       It provides the code to call all solvers from the Trilinos AztecOO package.

    """

    def __init__(self, tolerance=1e-10, criterion="default",
                 iterations=1000, precon=JacobiPreconditioner()):
        """
        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'unscaled', 'RHS', 'matrix', 'initial', 'solution', 'legacy'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.trilinos.preconditioners.preconditioner.Preconditioner
        """
        if self.__class__ is TrilinosAztecOOSolver:
            raise NotImplementedError("can't instantiate abstract base class")

        super(TrilinosAztecOOSolver, self).__init__(tolerance=tolerance, criterion=criterion,
                                                    iterations=iterations, precon=precon)

    def _adaptLegacyTolerance(self, L, x, b):
        return self._adaptInitialTolerance(L, x, b)

    def _adaptUnscaledTolerance(self, L, x, b):
        return (1., AztecOO.AZ_noscaled)

    def _adaptRHSTolerance(self, L, x, b):
        return (1., AztecOO.AZ_rhs)

    def _adaptMatrixTolerance(self, L, x, b):
        return (1., AztecOO.AZ_Anorm)

    def _adaptInitialTolerance(self, L, x, b):
        return (1., AztecOO.AZ_r0)

    def _adaptSolutionTolerance(self, L, x, b):
        return (1., AztecOO.AZ_sol)

    def _solve_(self, L, x, b):

        Solver = AztecOO.AztecOO(L, x, b)
        Solver.SetAztecOption(AztecOO.AZ_solver, self.solver)

##        Solver.SetAztecOption(AztecOO.AZ_kspace, 30)

        Solver.SetAztecOption(AztecOO.AZ_output, AztecOO.AZ_none)

        tolerance_scale, suite_criterion = self._adaptTolerance(L, x, b)

        rtol = self.scale_tolerance(self.tolerance, tolerance_scale)

        Solver.SetAztecOption(AztecOO.AZ_conv, suite_criterion)

        if self.preconditioner is None:
            Solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_none)
        else:
            self.preconditioner._applyToSolver(solver=Solver, matrix=L)

        self._log.debug("BEGIN solve")

        output = Solver.Iterate(self.iterations, rtol)

        self._log.debug("END solve")

        if self.preconditioner is not None:
            if hasattr(self.preconditioner, 'Prec'):
                del self.preconditioner.Prec

        status = Solver.GetAztecStatus()

        self._setConvergence(suite="trilinos",
                             code=int(status[AztecOO.AZ_why]),
                             iterations=int(status[AztecOO.AZ_its]),
                             tolerance_scale=tolerance_scale,
                             residual=status[AztecOO.AZ_r],
                             scaled_residual=status[AztecOO.AZ_scaled_r],
                             convergence_residual=status[AztecOO.AZ_rec_r],
                             solve_time=status[AztecOO.AZ_solve_time],
                             Aztec_version=status[AztecOO.AZ_Aztec_version])

        self.convergence.warn()

        return output
