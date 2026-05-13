__docformat__ = 'restructuredtext'

from scipy.sparse.linalg import cgs

from fipy.solvers.scipy.scipyKrylovSolver import ScipyKrylovSolver

__all__ = ["LinearCGSSolver"]

class LinearCGSSolver(ScipyKrylovSolver):
    """Interface to the conjugate gradient squared (:term:`CGS`) solver in
    :ref:`SciPy`.

    No preconditioning by default.
    """

    solveFnc = staticmethod(cgs)
