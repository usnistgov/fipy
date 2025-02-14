from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from petsc4py import PETSc

from fipy.tools.timer import Timer
from .petscSolver import PETScSolver
from .preconditioners.defaultPreconditioner import DefaultPreconditioner

from fipy.tools import numerix

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
                 divergence_tolerance=None,
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

    def _convergenceTest(self, ksp, its, rnorm, scale):
        """Replace `KSPConvergedDefault()` with custom normalization

        Modeled on `KSPConvergedDefault()
        <https://gitlab.com/petsc/petsc/-/blob/main/src/ksp/ksp/interface/iterativ.c#L1512>`_.
        Simplisitically, converged if `rnorm <= rtol * scale`.

        It would be much nicer (and less expensive) if they would just let
        you specify how to calculate rnorm0!
        `KSPConvergedDefaultSetUIRNorm()` isn't exposed to petsc4py, and it
        wouldn't help with "unscaled" or "matrix".

        Parameters
        ----------
        ksp : ~petsc4py.PETSc.KSP
            Krylov solver object.
        its : int
            Number of iterations performed.
        rnorm : float
            Norm of the residual to test.
        scale : float
            How to interpret magnitude of tolerance.

        Returns
        -------
        bool
            Whether solution is converged.
        """

        reason = PETSc.KSP.ConvergedReason.CONVERGED_ITERATING

        min_it = 0

        if numerix.isnan(rnorm) or numerix.isinf(rnorm):
            pcreason = PCReduceFailedReason
            if pcreason:
                reason = PETSc.KSP.ConvergedReason.DIVERGED_PCSETUP_FAILED
                self._log.debug("Linear solver pcsetup fails, "
                                "declaring divergence")
            else:
                reason = PETSc.KSP.ConvergedReason.DIVERGED_NANORINF
                self._log.debug("Linear solver has created a not a number (NaN) "
                                "as the residual norm, declaring divergence")
        elif its >= min_it:
            if rnorm <= max(ksp.rtol * scale, ksp.atol):
                if rnorm < ksp.atol:
                    reason = PETSc.KSP.ConvergedReason.CONVERGED_ATOL
                    self._log.debug("Linear solver has converged. "
                                    "Residual norm {rnorm:14.12e} is less "
                                    "than absolute tolerance {atol:14.12e} "
                                    "at iteration {its}".format(rnorm=rnorm,
                                                                atol=ksp.atol,
                                                                its=its))
                else:
                    reason = PETSc.KSP.ConvergedReason.CONVERGED_RTOL
                    self._log.debug("Linear solver has converged. "
                                    "Residual norm {rnorm:14.12e} is less "
                                    "than relative tolerance {rtol:14.12e} "
                                    "times residual scale {scale:14.12e} "
                                    "at iteration {its}".format(rnorm=rnorm,
                                                                rtol=ksp.rtol,
                                                                scale=scale,
                                                                its=its))
        elif rnorm >= ksp.dtol * scale:
            reason = PETSc.KSP.ConvergedReason.DIVERGED_DTOL
            self._log.debug("Linear solver is diverging. "
                            "Residual scale {scale:14.12e}, "
                            "current residual norm {rnorm:14.12e} "
                            "at iteration {its}".format(rnorm=rnorm,
                                                        scale=scale,
                                                        its=its))

        return reason

    def _adaptUnscaledTolerance(self, L, x, b):
        return (1., self._convergenceTest)

    def _adaptRHSTolerance(self, L, x, b):
        return (1., PETSc.KSP.NormType.UNPRECONDITIONED)

    def _adaptMatrixTolerance(self, L, x, b):
        return (self._matrixNorm(L, x, b), self._convergenceTest)

    def _adaptInitialTolerance(self, L, x, b):
        return (self._residualNorm(L, x, b), self._convergenceTest)

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
        if suite_criterion == self._convergenceTest:
            ksp.setConvergenceTest(suite_criterion,
                                   kargs=dict(scale=tolerance_scale))
            suite_criterion = PETSc.KSP.NormType.PRECONDITIONED
            rtol = self.tolerance
            divtol = self.divergence_tolerance
        else:
            rtol, divtol = (self.scale_tolerance(tol, tolerance_scale)
                            for tol in (self.tolerance,
                                        self.divergence_tolerance))

        ksp.setTolerances(rtol=rtol,
                          atol=self.absolute_tolerance,
                          divtol=divtol,
                          max_it=self.iterations)
        ksp.setNormType(suite_criterion)

        self._log.debug("BEGIN precondition")

        with Timer() as t:
            if self.preconditioner is None:
                ksp.getPC().setType("none")
            else:
                self.preconditioner._applyToSolver(solver=ksp, matrix=L)

        self._log.debug("END precondition - {} ns".format(t.elapsed))

        ksp.setFromOptions()

        self._log.debug("BEGIN solve")

        with Timer() as t:
            ksp.solve(b, x)

        self._log.debug("END solve - {} ns".format(t.elapsed))

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
