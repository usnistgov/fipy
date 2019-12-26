from __future__ import division
from builtins import range
from past.utils import old_div
__docformat__ = 'restructuredtext'

import os

from petsc4py import PETSc

from fipy.solvers.petsc.petscSolver import PETScSolver

__all__ = ["LinearLUSolver"]

class LinearLUSolver(PETScSolver):

    """
    The `LinearLUSolver` is an interface to the LU preconditioner in PETSc.
    A direct solve is performed.

    """
      
    def __init__(self, tolerance=1e-10, iterations=10, precon="lu"):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: *Ignored*. 

        """
        PETScSolver.__init__(self, tolerance=tolerance,
                             iterations=iterations, precon="lu")

    def _solve_(self, L, x, b):
        ksp = PETSc.KSP()
        ksp.create(PETSc.COMM_WORLD)
        ksp.setType("preonly")
        ksp.getPC().setType(self.preconditioner)
        # TODO: SuperLU invoked with PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU)
        #       see: http://www.mcs.anl.gov/petsc/petsc-dev/src/ksp/ksp/examples/tutorials/ex52.c.html
        # PETSc.PC().setFactorSolverType("superlu")
        
        L.assemblyBegin()
        L.assemblyEnd()
        ksp.setOperators(L)
        ksp.setFromOptions()
        
        for iteration in range(self.iterations):
            errorVector = L * x - b
            tol = errorVector.norm()
            
            if iteration == 0:
                tol0 = tol
                
            if (tol / tol0) <= self.tolerance:
                break
                
            xError = x.copy()

            ksp.solve(errorVector, xError)
            x -= xError
            
        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            from fipy.tools.debug import PRINT
#             L.view()
#             b.view()
            PRINT('solver:', ksp.type)
            PRINT('precon:', ksp.getPC().type)
            PRINT('iterations: %d / %d' % (iteration+1, self.iterations))
            PRINT('residual:', errorVector.norm(1))

