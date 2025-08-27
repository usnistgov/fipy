__docformat__ = 'restructuredtext'

from fipy.solvers.petsc.petscKrylovSolver import PETScKrylovSolver

__all__ = ["LinearCGSolver", "LinearPCGSolver"]

class LinearCGSolver(PETScKrylovSolver):

    """Interface to the conjugate gradient (:term:`CG`) solver in
    :ref:`PETSC`.
    """

    solver = 'cg'

    def _canSolveAsymmetric(self):
        return False

LinearPCGSolver = LinearCGSolver
