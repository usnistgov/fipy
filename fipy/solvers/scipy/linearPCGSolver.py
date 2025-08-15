from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from scipy.sparse.linalg import cg

from fipy.solvers.scipy.scipyKrylovSolver import ScipyKrylovSolver

__all__ = ["LinearPCGSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearPCGSolver(ScipyKrylovSolver):
    """Interface to the preconditioned conjugate gradient (:term:`PCG`) solver
    in :ref:`Scipy`.

    No preconditioning by default.
    """

    solveFnc = staticmethod(cg)

    def _canSolveAsymmetric(self):
        return False
