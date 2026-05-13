__docformat__ = 'restructuredtext'

from scipy.sparse.linalg import cg

from fipy.solvers.scipy.scipyKrylovSolver import ScipyKrylovSolver

__all__ = ["LinearCGSolver", "LinearPCGSolver"]

class LinearCGSolver(ScipyKrylovSolver):
    """Interface to the conjugate gradient (:term:`CG`) solver
    in :ref:`Scipy`.

    No preconditioning by default.
    """

    solveFnc = staticmethod(cg)

    def _canSolveAsymmetric(self):
        return False

LinearPCGSolver = LinearCGSolver
