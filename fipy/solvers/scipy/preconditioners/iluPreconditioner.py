# Adapted from https://stackoverflow.com/q/46876951/2019542

from __future__ import unicode_literals
from builtins import object
__docformat__ = 'restructuredtext'

__all__ = ["ILUPreconditioner"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

from scipy.sparse.linalg import LinearOperator, spilu

from .preconditioner import Preconditioner

class ILUPreconditioner(Preconditioner):
    
    def _applyToMatrix(self, A):
        ilu = spilu(A.tocsc())
        Mx = lambda x: ilu.solve(x)

        return LinearOperator(A.shape, Mx)
