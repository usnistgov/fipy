from __future__ import division
from builtins import range
from past.utils import old_div
__docformat__ = 'restructuredtext'

from petsc4py import PETSc

from .petscSolver import PETScSolver
from .preconditioners.luPreconditioner import LUPreconditioner

__all__ = ["LinearLUSolver"]

class LinearLUSolver(PETScSolver):

    """
    The `LinearLUSolver` is an interface to the LU preconditioner in PETSc.
    A direct solve is performed.

    """

    def __init__(self, tolerance="default", criterion="default",
                 iterations=10, precon=None):
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
        precon
            *ignored*
        """
        super(LinearLUSolver, self).__init__(tolerance=tolerance,
                                             criterion=criterion,
                                             iterations=iterations,
                                             precon=LUPreconditioner())

    def _adaptLegacyTolerance(self, L, x, b):
        return self._adaptInitialTolerance(L, x, b)

    def _adaptUnscaledTolerance(self, L, x, b):
        return (1., None)

    def _adaptRHSTolerance(self, L, x, b):
        return (self._rhsNorm(L, x, b), None)

    def _adaptMatrixTolerance(self, L, x, b):
        return (self._matrixNorm(L, x, b), None)

    def _adaptInitialTolerance(self, L, x, b):
        return (self._residualNorm(L, x, b), None)

    def _solve_(self, L, x, b):
        ksp = PETSc.KSP()
        ksp.create(PETSc.COMM_WORLD)
        ksp.setType("preonly")
        self.preconditioner._applyToSolver(solver=ksp, matrix=L)
        # TODO: SuperLU invoked with PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU)
        #       see: http://www.mcs.anl.gov/petsc/petsc-dev/src/ksp/ksp/examples/tutorials/ex52.c.html
        # PETSc.PC().setFactorSolverType("superlu")

        L.assemble()
        ksp.setOperators(L)
        ksp.setFromOptions()

        tolerance_scale, _ = self._adaptTolerance(L, x, b)

        self._log.debug("BEGIN solve")

        for iteration in range(self.iterations):
            residualVector, residual = self._residualVectorAndNorm(L, x, b)

            if residual <= self.tolerance * tolerance_scale:
                break

            xError = x.copy()

            ksp.solve(residualVector, xError)

            x -= xError

        self._log.debug("END solve")

        self._setConvergence(suite="petsc",
                             code=PETSc.KSP.ConvergedReason.CONVERGED_ITS,
                             iterations=iteration+1,
                             residual=residual,
                             ksp_solver=ksp.type,
                             ksp_precon=ksp.getPC().type)

        self.convergence.warn()
