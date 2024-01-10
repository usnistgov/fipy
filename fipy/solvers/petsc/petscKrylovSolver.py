from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from petsc4py import PETSc

from .petscSolver import PETScSolver
from .preconditioners.defaultPreconditioner import DefaultPreconditioner

__all__ = ["PETScKrylovSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class PETScKrylovSolver(PETScSolver):

    """
    .. attention:: This class is abstract, always create one of its subclasses.
       It provides the code to call all Krylov solvers from the PETSc package.

    """

    DEFAULT_PRECONDITIONER = DefaultPreconditioner

    def __init__(self, tolerance="default",
                 absolute_tolerance=None,
                 divergence_tolerance=1e100,
                 criterion="default",
                 iterations="default", precon="default"):
        """
        Parameters
        ----------
        tolerance : float
            Required relative error tolerance.
        absolute_tolerance : float
            Required absolute error tolerance.
        divergence_tolerance : float
            Required divergence error tolerance.
        criterion : {'default', 'unscaled', 'RHS', 'matrix', 'initial', 'preconditioned', 'natural', 'legacy'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : ~fipy.solvers.petsc.preconditioners.petscPreconditioner.PETScPreconditioner, optional
            Preconditioner to apply to the matrix.  A value of None means
            to perform an unpreconditioned solve.  (default:
            :class:`~fipy.solvers.petsc.preconditioners.defaultPreconditioner.DefaultPreconditioner`).
        """
        if self.__class__ is PETScKrylovSolver:
            raise NotImplementedError("can't instantiate abstract base class")

        self.absolute_tolerance = absolute_tolerance
        self.divergence_tolerance = divergence_tolerance
        super(PETScKrylovSolver, self).__init__(tolerance=tolerance, criterion=criterion,
                                                iterations=iterations, precon=precon)

    def _adaptLegacyTolerance(self, L, x, b):
        return (1., PETSc.KSP.NormType.DEFAULT)

    def _adaptUnscaledTolerance(self, L, x, b):
        factor = 1. / self._rhsNorm(L, x, b)
        return (factor, PETSc.KSP.NormType.UNPRECONDITIONED)

    def _adaptRHSTolerance(self, L, x, b):
        return (1., PETSc.KSP.NormType.UNPRECONDITIONED)

    def _adaptMatrixTolerance(self, L, x, b):
        factor = self._matrixNorm(L, x, b) / self._rhsNorm(L, x, b)
        return (factor, PETSc.KSP.NormType.UNPRECONDITIONED)

    def _adaptInitialTolerance(self, L, x, b):
        factor = self._residualNorm(L, x, b) / self._rhsNorm(L, x, b)
        return (factor, PETSc.KSP.NormType.UNPRECONDITIONED)

    def _adaptPreconditionedTolerance(self, L, x, b):
        return (1., PETSc.KSP.NormType.PRECONDITIONED)

    def _adaptNaturalTolerance(self, L, x, b):
        return (1., PETSc.KSP.NormType.NATURAL)

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
        ksp.create(L.comm)
        ksp.setType(self.solver)
        if self.criterion != "legacy":
            ksp.setInitialGuessNonzero(True)

        L.assemble()
        ksp.setOperators(L)

        tolerance_scale, suite_criterion = self._adaptTolerance(L, x, b)

        rtol, divtol = (self.scale_tolerance(tol, tolerance_scale)
                        for tol in (self.tolerance,
                                    self.divergence_tolerance))

        ksp.setTolerances(rtol=rtol,
                          atol=self.absolute_tolerance,
                          divtol=divtol,
                          max_it=self.iterations)
        ksp.setNormType(suite_criterion)

        self._log.debug("BEGIN precondition")

        if self.preconditioner is None:
            ksp.getPC().setType("none")
        else:
            self.preconditioner._applyToSolver(solver=ksp, matrix=L)

        self._log.debug("END precondition")

        ksp.setFromOptions()

        self._log.debug("BEGIN solve")

        ksp.solve(b, x)

        self._log.debug("END solve")

        self._setConvergence(suite="petsc",
                             code=ksp.reason,
                             iterations=ksp.its,
                             tolerance_scale=tolerance_scale,
                             # "The residual value that is tested may be an approximation"
                             # https://petsc.org/release/petsc4py/reference/petsc4py.PETSc.KSP.html#petsc4py.PETSc.KSP.setConvergenceTest
                             residual=ksp.norm,
                             ksp_solver=ksp.type,
                             ksp_precon=ksp.getPC().type,
                             ksp_norm_type=ksp.norm_type)

        ksp.destroy()

        self.convergence.warn()

        return x
