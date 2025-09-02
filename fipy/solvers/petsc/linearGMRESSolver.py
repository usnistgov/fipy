__docformat__ = 'restructuredtext'

from fipy.solvers.petsc.petscKrylovSolver import PETScKrylovSolver

__all__ = ["LinearGMRESSolver"]

class LinearGMRESSolver(PETScKrylovSolver):

    """Interface to the generalized minimal residual (:term:`GMRES`) solver
    in :ref:`PETSC`.
    """
    
    solver = 'gmres'
