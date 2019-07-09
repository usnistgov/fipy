__docformat__ = 'restructuredtext'

from fipy.solvers.petsc.petscKrylovSolver import PETScKrylovSolver

__all__ = ["LinearBicgSolver"]

class LinearBicgSolver(PETScKrylovSolver):

    """
    The `LinearBicgSolver` is an interface to the biconjugate gradient solver in
    PETSc, using no preconditioner by default.
    """
      
    solver = 'bicg'
