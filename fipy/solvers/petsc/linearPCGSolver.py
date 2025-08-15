__docformat__ = 'restructuredtext'

from fipy.solvers.petsc.petscKrylovSolver import PETScKrylovSolver

__all__ = ["LinearPCGSolver"]

class LinearPCGSolver(PETScKrylovSolver):

    """Interface to the conjugate graident (:term:`CG`) solver in
    :ref:`PETSC`.
    """
      
    solver = 'cg'
