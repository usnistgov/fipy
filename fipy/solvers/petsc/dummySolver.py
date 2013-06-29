__docformat__ = 'restructuredtext'

from fipy.solvers.petsc.petscSolver import PETScSolver

__all__ = ["DummySolver"]

class DummySolver(PETScSolver):

    """Solver that doesn't do anything.
    
    PETSc is intolerant of having zeros on the diagonal
    """
      
    def _solve_(self, L, x, b):
        pass
