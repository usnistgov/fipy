__docformat__ = 'restructuredtext'

import os

from petsc4py import PETSc

from fipy.solvers.petsc.petscSolver import PETScSolver

__all__ = ["PETScKrylovSolver"]

_reason = {1: "KSP_CONVERGED_RTOL_NORMAL",
           9: "KSP_CONVERGED_ATOL_NORMAL",
           2: "KSP_CONVERGED_RTOL",
           3: "KSP_CONVERGED_ATOL",
           4: "KSP_CONVERGED_ITS",
           5: "KSP_CONVERGED_CG_NEG_CURVE",
           6: "KSP_CONVERGED_CG_CONSTRAINED",
           7: "KSP_CONVERGED_STEP_LENGTH",
           8: "KSP_CONVERGED_HAPPY_BREAKDOWN",
           -2: "KSP_DIVERGED_NULL",
           -3: "KSP_DIVERGED_ITS",
           -4: "KSP_DIVERGED_DTOL",
           -5: "KSP_DIVERGED_BREAKDOWN",
           -6: "KSP_DIVERGED_BREAKDOWN_BICG",
           -7: "KSP_DIVERGED_NONSYMMETRIC",
           -8: "KSP_DIVERGED_INDEFINITE_PC",
           -9: "KSP_DIVERGED_NANORINF",
           -10: "KSP_DIVERGED_INDEFINITE_MAT",
           -11: "KSP_DIVERGED_PC_FAILED",
           0: "KSP_CONVERGED_ITERATING"}

class PETScKrylovSolver(PETScSolver):

    """
    .. attention:: This class is abstract, always create one of its subclasses.
       It provides the code to call all Krylov solvers from the PETSc package.

    """
      
    def __init__(self, tolerance=1e-10, iterations=1000, precon=None):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner to use (string). 

        """
        if self.__class__ is PETScKrylovSolver:
            raise NotImplementedError("can't instantiate abstract base class")
            
        PETScSolver.__init__(self, tolerance=tolerance,
                             iterations=iterations, precon=precon)

    def _solve_(self, L, x, b):
        ksp = PETSc.KSP()
        ksp.create(L.comm)
        ksp.setType(self.solver)
        if self.preconditioner is not None:
            ksp.getPC().setType(self.preconditioner)
        ksp.setTolerances(rtol=self.tolerance, max_it=self.iterations)
        L.assemble()
        ksp.setOperators(L)
        ksp.setFromOptions()
        ksp.solve(b, x)

        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            from fipy.tools.debug import PRINT
#             L.view()
#             b.view()
            PRINT('solver:', ksp.type)
            PRINT('precon:', ksp.getPC().type)
            PRINT('convergence: %s' % _reason[ksp.reason])
            PRINT('iterations: %d / %d' % (ksp.its, self.iterations))
            PRINT('norm:', ksp.norm)
            PRINT('norm_type:', ksp.norm_type)
