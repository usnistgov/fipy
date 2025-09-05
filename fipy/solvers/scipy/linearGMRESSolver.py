__docformat__ = 'restructuredtext'

from scipy.sparse.linalg import gmres

from fipy.solvers.scipy.scipyKrylovSolver import ScipyKrylovSolver

__all__ = ["LinearGMRESSolver"]

class LinearGMRESSolver(ScipyKrylovSolver):
    """Interface to the Generalized Minimum RESidual (:term:`GMRES`) solver
    in :ref:`SciPy`.

    No preconditioning by default.
    """

    solveFnc = staticmethod(gmres)

    def _doSolve(self, *args, **kwargs):
        return self.solveFnc(*args, **kwargs, callback_type='legacy')
