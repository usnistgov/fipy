from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from scipy.sparse.linalg import gmres

from fipy.solvers.scipy.scipyKrylovSolver import ScipyKrylovSolver

__all__ = ["LinearGMRESSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearGMRESSolver(ScipyKrylovSolver):
    """Interface to the Generalized Minimum RESidual (:term:`GMRES`) solver
    in :ref:`SciPy`.

    No preconditioning by default.
    """

    solveFnc = staticmethod(gmres)

    def _doSolve(self, *args, **kwargs):
        return self.solveFnc(*args, **kwargs, callback_type='legacy')
