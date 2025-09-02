from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.solvers.scipy.scipyKrylovSolver import ScipyKrylovSolver
from scipy.sparse.linalg import bicgstab

__all__ = ["LinearBicgstabSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LinearBicgstabSolver(ScipyKrylovSolver):
    """Interface to the Biconjugate Gradient (Stabilized) (:term:`BiCGSTAB`)
    solver in :ref:`Scipy`.

    No preconditioning by default.
    """

    solveFnc = staticmethod(bicgstab)
