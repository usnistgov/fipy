__docformat__ = 'restructuredtext'

import os

from petsc4py import PETSc

from fipy.solvers.petsc.petscSolver import PETScSolver

__all__ = ["PETScKrylovSolver"]

class PETScKrylovSolver(PETScSolver):

    """
    .. attention:: This class is abstract, always create on of its subclasses. It provides the code to call all solvers from the Trilinos AztecOO package.

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
        L.assemblyBegin()
        L.assemblyEnd()
        ksp.setOperators(L)
        ksp.setFromOptions()
        ksp.solve(b, x)
