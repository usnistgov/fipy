__docformat__ = 'restructuredtext'

from fipy.solvers.petsc.petscKrylovSolver import PETScKrylovSolver

__all__ = ["LinearCGSSolver"]

class LinearCGSSolver(PETScKrylovSolver):

    """
    The `LinearCGSSolver` is an interface to the conjugate gradient squared
    solver in PETSc, using no preconditioner by default.
    """
      
    solver = 'cgs'
