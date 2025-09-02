__docformat__ = 'restructuredtext'

from fipy.solvers.petsc.petscKrylovSolver import PETScKrylovSolver

__all__ = ["LinearCGSSolver"]

class LinearCGSSolver(PETScKrylovSolver):

    """Interface to the conjugate gradient squared (:term:`CGS`) solver in
    :ref:`PETSc`.
    """
      
    solver = 'cgs'
