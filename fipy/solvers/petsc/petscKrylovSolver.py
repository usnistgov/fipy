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

    criteria = {
        "default": PETSc.KSP.NormType.DEFAULT
        "preconditioned": PETSc.KSP.NormType.PRECONDITIONED,
        "initial": PETSc.KSP.NormType.UNPRECONDITIONED,
        "natural": PETSc.KSP.NormType.NATURAL,
    }

    def __init__(self, tolerance=1e-10, criterion="default",
                 iterations=1000, precon=None):
        """
        Parameters
        ----------
        tolerance : float
            Required error tolerance.
        criterion : {'default', 'initial', 'preconditioned', 'natural'}
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

    def _solve_(self, L, x, b):
        ksp = PETSc.KSP()
        ksp.create(L.comm)
        ksp.setType(self.solver)
        if self.preconditioner is not None:
            ksp.getPC().setType(self.preconditioner)
        ksp.setTolerances(rtol=self.tolerance, max_it=self.iterations)
        ksp.setNormType(self.criteria[self.criterion])
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
