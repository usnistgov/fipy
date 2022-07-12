from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from petsc4py import PETSc

from .petscSolver import PETScSolver

__all__ = ["PETScKrylovSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class PETScKrylovSolver(PETScSolver):

    """
    .. attention:: This class is abstract, always create one of its subclasses.
       It provides the code to call all Krylov solvers from the PETSc package.

    """

    def __init__(self, tolerance=1e-10, criterion="default",
                 iterations=1000, precon=None):
        """
        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'unscaled', 'RHS', 'matrix', 'initial', 'preconditioned', 'natural'}
            Interpretation of ``tolerance``.
            See :ref:`CONVERGENCE` for more information.
        iterations : int
            Maximum number of iterative steps to perform.
        precon : str
            Preconditioner to use.

        """
        if self.__class__ is PETScKrylovSolver:
            raise NotImplementedError("can't instantiate abstract base class")
            
        PETScSolver.__init__(self, tolerance=tolerance, criterion=criterion,
                             iterations=iterations, precon=precon)

    def _adaptDefaultTolerance(self, L, x, b):
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
        ksp = PETSc.KSP()
        ksp.create(L.comm)
        ksp.setType(self.solver)
        if self.preconditioner is None:
            ksp.getPC().setType('none')
        else:
            self.preconditioner._applyToSolver(solver=ksp, matrix=L)

        tolerance_factor, suite_criterion = self._adaptTolerance(L, x, b)
        ksp.setTolerances(rtol=self.tolerance * tolerance_factor, max_it=self.iterations)
        ksp.setNormType(suite_criterion)

        L.assemble()
        ksp.setOperators(L)
        ksp.setFromOptions()
        ksp.solve(b, x)

        self._setConvergence(suite="petsc",
                             code=ksp.reason,
                             iterations=ksp.its,
                             residual=ksp.norm,
                             ksp_solver=ksp.type,
                             ksp_precon=ksp.getPC().type,
                             ksp_norm_type=ksp.norm_type)

        ksp.destroy()

        self.convergence.warn()
