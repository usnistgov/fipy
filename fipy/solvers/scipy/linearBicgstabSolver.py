__docformat__ = 'restructuredtext'

from fipy.solvers.scipy.scipyKrylovSolver import ScipyKrylovSolver
from scipy.sparse.linalg import bicgstab

__all__ = ["LinearBicgstabSolver"]

class LinearBicgstabSolver(ScipyKrylovSolver):
    """Interface to the Biconjugate Gradient (Stabilized) (:term:`BiCGSTAB`)
    solver in :ref:`Scipy`.

    No preconditioning by default.
    """

    solveFnc = staticmethod(bicgstab)
