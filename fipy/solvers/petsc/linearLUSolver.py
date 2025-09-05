__docformat__ = 'restructuredtext'

from petsc4py import PETSc

from fipy.tools.timer import Timer
from .petscSolver import PETScSolver
from .preconditioners.luPreconditioner import LUPreconditioner

__all__ = ["LinearLUSolver"]

class LinearLUSolver(PETScSolver):

    """Interface to the :term:`LU` preconditioner in :ref:`PETSC`.

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
        """Solve system of equations posed for PETSc

        Parameters
        ----------
        L : PETSc.Mat
            Sparse matrix
        x : PETSc.Vec
            Solution variable as ghosted vector
        b : PETSc.Vec
            Right-hand side as ghosted vector

        Returns
        -------
        x : PETSc.Vec
            Solution variable as ghosted vector
        """
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

        with Timer() as t:
            for iteration in range(self.iterations):
                residualVector, residual = self._residualVectorAndNorm(L, x, b)

                if residual <= self.tolerance * tolerance_scale:
                    break

                xError = x.copy()

                ksp.solve(residualVector, xError)

                x -= xError

        self._log.debug("END solve - {} ns".format(t.elapsed))

        self._setConvergence(suite="petsc",
                             code=PETSc.KSP.ConvergedReason.CONVERGED_ITS,
                             iterations=iteration+1,
                             residual=residual,
                             ksp_solver=ksp.type,
                             ksp_precon=ksp.getPC().type)

        self.convergence.warn()

        return x
