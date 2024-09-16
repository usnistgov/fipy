from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = ["ILUPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

from scipy.sparse.linalg import LinearOperator, spilu

from .scipyPreconditioner import ScipyPreconditioner

class ILUPreconditioner(ScipyPreconditioner):
    """ILU preconditioner for :class:`~fipy.solvers.scipy.scipySolver.ScipySolver`.

    Wrapper class for :func:`scipy.sparse.linalg.spilu`.
    Adapted from https://stackoverflow.com/q/46876951/2019542.
    """
    
    def _applyToMatrix(self, matrix):
        ilu = spilu(matrix.tocsc())
        Mx = lambda x: ilu.solve(x)

        return LinearOperator(matrix.shape, Mx), matrix
