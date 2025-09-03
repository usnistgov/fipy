__docformat__ = 'restructuredtext'

from PyTrilinos import AztecOO

from .trilinosSolver import TrilinosSolver
from .preconditioners.jacobiPreconditioner import JacobiPreconditioner
from fipy.tools.timer import Timer

__all__ = ["TrilinosAztecOOSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class TrilinosAztecOOSolver(TrilinosSolver):

    """
    .. attention:: This class is abstract, always create on of its subclasses.
       It provides the code to call all solvers from the Trilinos AztecOO package.

    """

    DEFAULT_PRECONDITIONER = JacobiPreconditioner

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
        """Solve system of equations posed for PyTrilinos

        Parameters
        ----------
        L : Epetra.CrsMatrix
            Sparse matrix
        x : Epetra.Vector
            Solution variable as non-ghosted vector
        b : Epetra.Vector
            Right-hand side as non-ghosted vector

        Returns
        -------
        x : Epetra.Vector
            Solution variable as non-ghosted vector
        """

        solver = AztecOO.AztecOO(L, x, b)
        solver.SetAztecOption(AztecOO.AZ_solver, self.solver)

##        solver.SetAztecOption(AztecOO.AZ_kspace, 30)

        solver.SetAztecOption(AztecOO.AZ_output, AztecOO.AZ_none)

        tolerance_scale, suite_criterion = self._adaptTolerance(L, x, b)

        rtol = self.scale_tolerance(self.tolerance, tolerance_scale)

        solver.SetAztecOption(AztecOO.AZ_conv, suite_criterion)

        self._log.debug("BEGIN precondition")

        with Timer() as t:
            if self.preconditioner is None:
                solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_none)
            else:
                self.preconditioner._applyToSolver(solver=solver, matrix=L)

        self._log.debug("END precondition - {} ns".format(t.elapsed))

        self._log.debug("BEGIN solve")

        with Timer() as t:
            solver.Iterate(self.iterations, rtol)

        self._log.debug("END solve - {} ns".format(t.elapsed))

        if self.preconditioner is not None:
            if hasattr(self.preconditioner, 'Prec'):
                del self.preconditioner.Prec

        status = solver.GetAztecStatus()

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

        return x
