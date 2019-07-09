__docformat__ = 'restructuredtext'

from fipy.solvers.petsc.petscKrylovSolver import PETScKrylovSolver

__all__ = ["LinearPCGSolver"]

class LinearPCGSolver(PETScKrylovSolver):

    """
    The `LinearPCGSolver` is an interface to the cg solver in PETSc,
    using no preconditioner by default.

    """
      
    solver = 'cg'
