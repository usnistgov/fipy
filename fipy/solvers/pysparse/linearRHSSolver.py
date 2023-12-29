from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from .pysparseSolver import PysparseSolver

class LinearRHSSolver(PysparseSolver):
    """Wrapper for solvers that normalize the residual by the right-hand-side.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def _legacyNorm(self, L, x, b):
        return self._rhsNorm(L, x, b)