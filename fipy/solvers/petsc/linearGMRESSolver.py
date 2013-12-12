__docformat__ = 'restructuredtext'

from fipy.solvers.petsc.petscKrylovSolver import PETScKrylovSolver

__all__ = ["LinearGMRESSolver"]

class LinearGMRESSolver(PETScKrylovSolver):

    """
    The `LinearGMRESSolver` is an interface to the GMRES solver in PETSc,
    using no preconditioner by default.

    """
    
    solver = 'gmres'
